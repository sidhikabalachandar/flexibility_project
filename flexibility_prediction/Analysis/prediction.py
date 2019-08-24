import pandas as pd
import time
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.tree import DecisionTreeRegressor
from sklearn.preprocessing import StandardScaler
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

dict_classifiers = {
    "Linear Regression": LinearRegression(),
    "Polynomial Regression": LinearRegression(),
#     "Neural Network": MLPRegressor(hidden_layer_sizes=(8,8,8), activation='relu', solver='adam', max_iter=500),
    "Decision Tree": DecisionTreeRegressor(random_state = 0),
    "Decision Tree with Max Depth = 03": DecisionTreeRegressor(max_depth=3, random_state = 0),
    "Decision Tree with Max Depth = 05": DecisionTreeRegressor(max_depth=5, random_state = 0),
    "Decision Tree with Max Depth = 10": DecisionTreeRegressor(max_depth=10, random_state = 0),
    "Decision Tree with Max Depth = 15": DecisionTreeRegressor(max_depth=15, random_state = 0),
#     "Random Forest": RandomForestRegressor(n_estimators=100, random_state=0)
}


def batch_classify(X_train, Y_train, X_test, Y_test, verbose=True, include_y_pred=False):
    """
    This method, takes as input the X, Y matrices of the Train and Test set.
    And fits them on all of the Classifiers specified in the dict_classifier.
    The trained models, and accuracies are saved in a dictionary. The reason to use a dictionary
    is because it is very easy to save the whole dictionary with the pickle module.

    Usually, the SVM, Random Forest and Gradient Boosting Classifier take quiet some time to train.
    So it is best to train them on a smaller dataset first and
    decide whether you want to comment them out or not based on the test accuracy score.
    """

    dict_models = {}
    for classifier_name, classifier in list(dict_classifiers.items()):
        coefs = None
        if classifier_name == "Polynomial Regression":
            poly = PolynomialFeatures(degree=2)
            saved_X_train = X_train
            saved_X_test = X_test
            X_train = poly.fit_transform(X_train)
            X_test = poly.fit_transform(X_test)
        if classifier_name == "Random Forest":
            sc = StandardScaler()
            X_train = sc.fit_transform(X_train)
            X_test = sc.transform(X_test)
        t_start = time.process_time()
        classifier.fit(X_train, Y_train)
        t_end = time.process_time()

        if classifier_name == "Linear Regression":
            coefs = classifier.coef_

        t_diff = t_end - t_start
        train_score = classifier.score(X_train, Y_train)
        Y_pred = classifier.predict(X_test)
        test_score = classifier.score(X_test, Y_test)

        if include_y_pred:
            dict_models[classifier_name] = {'model': classifier_name, 'train_score': train_score,
                                            'test_score': test_score, 'train_time': t_diff, 'y_pred': Y_pred}
        else:
            dict_models[classifier_name] = {'model': classifier_name, 'train_score': train_score,
                                            'test_score': test_score, 'train_time': t_diff, 'coef': coefs}
        if verbose:
            print("trained {c} in {f:.2f} s".format(c=classifier_name, f=t_diff))
        if classifier_name == "Polynomial Regression":
            X_train = saved_X_train
            X_test = saved_X_test
    return dict_models


def scores_calculator(Y_test, classifiers, model):
    score_dict = {}
    for classifier in classifiers:
        Y_pred_bin = model[classifier]['y_pred'] > 2
        Y_test_bin = Y_test > 2
        true_pos_count = 0
        true_neg_count = 0
        false_pos_count = 0
        false_neg_count = 0
        for i in range(len(Y_pred_bin)):
            if Y_pred_bin[i] == Y_test_bin[i]:
                if Y_pred_bin[i] == True:
                    true_pos_count += 1
                else:
                    true_neg_count += 1
            else:
                if Y_pred_bin[i] == True:
                    false_pos_count += 1
                else:
                    false_neg_count += 1

        accuracy = (true_pos_count + true_neg_count) / len(Y_pred_bin)
        if (true_pos_count + false_neg_count) != 0:
            recall = true_pos_count / (true_pos_count + false_neg_count)
        else:
            recall = None
        if (true_pos_count + false_pos_count) != 0:
            precision = true_pos_count / (true_pos_count + false_pos_count)
        else:
            precision = None

        score_dict[classifier] = {'Accuracy': accuracy,
                                  'Recall': recall,
                                  'Precision': precision}

    return score_dict

def load_data():
    files = os.listdir('../Data/rmsds')
    data = pd.read_csv("../Data/rmsds/5HT2B_rmsds.csv")
    for file in files:
        if file[-4:] == '.csv' and file != '5HT2B_rmsds.csv':
            fileData = pd.read_csv("../Data/rmsds/" + file)
            data = pd.concat([data, fileData], sort=False)
    data = data[data['secondary structure'] != -1]
    data = data[data['complete rmsd'] < 30]
    prots = sorted(data['protein'].unique())
    return (data, prots)

def run_pred(train_data, test_data, ignore_cols):
    X_train = train_data.drop(ignore_cols, axis=1).values
    Y_train = train_data['complete rmsd'].values

    X_test = test_data.drop(ignore_cols, axis=1).values
    Y_test = test_data['complete rmsd'].values

    X_true = X_train[Y_train > 2]
    Y_true = Y_train[Y_train > 2]

    X_train_expanded_50 = np.append(X_train, X_true, axis=0)
    Y_train_expanded_50 = np.append(Y_train, Y_true, axis=0)

    for i in range(4):
        X_train_expanded_50 = np.append(X_train_expanded_50, X_true, axis=0)
        Y_train_expanded_50 = np.append(Y_train_expanded_50, Y_true, axis=0)

    dict_models = batch_classify(X_train_expanded_50, Y_train_expanded_50, X_test, Y_test, verbose=False, include_y_pred=True)
    scores = scores_calculator(Y_test, list(dict_classifiers.keys()), dict_models)
    return (dict_models, scores)

def results(test_scores, prots):
    avg_scores = {}
    counts = {}
    tree_depth_05 = {}
    for classifier in dict_classifiers:
        avg_scores[classifier] = {'Accuracy': 0, 'Recall': 0, 'Precision': 0}
        counts[classifier] = {'Accuracy': 0, 'Recall': 0, 'Precision': 0}

    for i, scores in enumerate(test_scores):
        for classifier in scores:
            if classifier == 'Decision Tree with Max Depth = 05':
                tree_depth_05[prots[i]] = {}
            for score_type in scores[classifier]:
                if classifier == 'Decision Tree with Max Depth = 05':
                    tree_depth_05[prots[i]][score_type] = scores[classifier][score_type]
                if scores[classifier][score_type] != None:
                    avg_scores[classifier][score_type] += scores[classifier][score_type]
                    counts[classifier][score_type] += 1

    for classifier in avg_scores:
        for score_type in avg_scores[classifier]:
            avg_scores[classifier][score_type] /= counts[classifier][score_type]

    print(pd.DataFrame.from_dict(avg_scores, orient='index'))
    print(pd.DataFrame.from_dict(tree_depth_05, orient='index'))

def actual_rmsd_grapher(data, test_prot, start, target):
    selected = data[(data['protein'] == test_prot) & (data['start ligand'] == start) & (data['target ligand'] == target)]
    rmsd = np.array(selected['backbone rmsd'])
    rmsd[rmsd > 8] = 8
    rmsd_formatted = np.expand_dims(rmsd, axis=0)
    sns.set_context("talk", font_scale=1.0)

    plt.imshow(rmsd_formatted, cmap="plasma", interpolation='none')

    fig = plt.gcf()
    fig.set_size_inches(10, 3)

    ax = plt.gca()
    labels = [r + '\n' + str(i) for r, i in zip(selected['name'], selected['num'])]
    plt.xticks(np.arange(0, len(selected)), labels)
    plt.yticks([], [])
    ax.tick_params(axis=u'both', which=u'both', length=0)
    ax.tick_params(axis='both', labelsize=5.5)

    plt.colorbar()
    plt.title('F10, start ligand = 2J94, target ligand = 1Z6E')
    plt.show()

def pred_rmsd_grapher(dict_models, model, data, test_prot, start, target):
    selected = data[(data['protein'] == test_prot) & (data['start ligand'] == start) & (data['target ligand'] == target)]
    bin_pred_rmsd_formatted = np.expand_dims(dict_models[model]['y_pred'], axis=0)
    sns.set_context("talk", font_scale=1.0)
    plt.imshow(bin_pred_rmsd_formatted, cmap="plasma", interpolation='none')
    fig = plt.gcf()
    fig.set_size_inches(10, 3)

    ax = plt.gca()
    labels = [r + '\n' + str(i) for r, i in zip(selected['name'], selected['num'])]
    plt.xticks(np.arange(0, len(selected)), labels)
    plt.yticks([], [])
    ax.tick_params(axis=u'both', which=u'both', length=0)
    ax.tick_params(axis='both', labelsize=5.5)

    plt.colorbar()
    plt.title(model)
    plt.show()