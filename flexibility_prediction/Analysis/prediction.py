import pandas as pd
import time
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
import sklearn.metrics as metrics
import sklearn.utils as utils
import os


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
        t_start = time.process_time()
        classifier.fit(X_train, Y_train)
        t_end = time.process_time()

        t_diff = t_end - t_start
        train_score = classifier.score(X_train, Y_train)
        Y_pred = classifier.predict(X_test)
        test_score = classifier.score(X_test, Y_test)
        recall = metrics.recall_score(Y_test, Y_pred)
        precision = metrics.precision_score(Y_test, Y_pred)

        if include_y_pred:
            dict_models[classifier_name] = {'model': classifier_name, 'train_score': train_score,
                                            'test_score': test_score, 'recall_score': recall,
                                            'precision_score': precision, 'train_time': t_diff, 'y_pred': Y_pred}
        else:
            dict_models[classifier_name] = {'model': classifier_name, 'train_score': train_score,
                                            'test_score': test_score, 'recall_score': recall,
                                            'precision_score': precision, 'train_time': t_diff}
        if verbose:
            print("trained {c} in {f:.2f} s".format(c=classifier_name, f=t_diff))
    return dict_models


if __name__ == '__main__':
    files = os.listdir('../Data/rmsds')
    first = True
    data = pd.read_csv("../Data/rmsds/5HT2B_rmsds.csv")

    for file in files:
        if file[-4:] == '.csv' and file != '5HT2B_rmsds.csv':
            fileData = pd.read_csv("../Data/rmsds/" + file)
            data = pd.concat([data, fileData])
    data = data[data['secondary structure'] != -1]
    prots = data['protein'].unique()

    utils.shuffle(prots, random_state=7)

    train_prots = prots[:len(prots) * 66 // 100]
    test_prots = prots[len(prots) * 66 // 100:]

    train_data = data[data['protein'].isin(train_prots)]
    test_data = data[data['protein'].isin(test_prots)]

    X_train = train_data.drop(['protein', 'start ligand', 'target ligand', 'rmsd', 'normalized bfactor', 'res name'],
                              axis=1).values
    Y_train = train_data['rmsd'].values > 2

    X_test = test_data.drop(['protein', 'start ligand', 'target ligand', 'rmsd', 'normalized bfactor', 'res name'],
                            axis=1).values
    Y_test = test_data['rmsd'].values > 2

    X_true = X_train[Y_train == True]
    Y_true = Y_train[Y_train == True]

    dict_classifiers = {
        "Logistic Regression": LogisticRegression(),
        #     "Nearest Neighbors": KNeighborsClassifier(),
        #     "Linear SVM": SVC(),
        #     "Gradient Boosting Classifier": GradientBoostingClassifier(n_estimators=1000),
        "Decision Tree": tree.DecisionTreeClassifier(),
        "Decision Tree with max depth = 03": tree.DecisionTreeClassifier(max_depth=3),
        "Decision Tree with max depth = 05": tree.DecisionTreeClassifier(max_depth=5),
        "Decision Tree with max depth = 10": tree.DecisionTreeClassifier(max_depth=10),
        "Decision Tree with max depth = 15": tree.DecisionTreeClassifier(max_depth=15),
        "Random Forest": RandomForestClassifier(n_estimators=100),
        #     "Neural Net": MLPClassifier(alpha = 1),
        "Naive Bayes": GaussianNB(),
        # "AdaBoost": AdaBoostClassifier(),
        # "QDA": QuadraticDiscriminantAnalysis(),
        # "Gaussian Process": GaussianProcessClassifier()
    }