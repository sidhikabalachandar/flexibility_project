def grouper(n, iterable, iterable_2):
    iterable = list(iterable)
    iterable_2 = list(iterable_2)
    out = []
    out_2 = []
    for i in range(0, len(iterable), n):
        out += [iterable[i: i+n]]
        out_2 += [iterable_2[i: i+n]]
    return (out, out_2)
