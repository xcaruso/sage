class PrettyPrint():
    def __init__(self, X, top=7, bottom=3):
        self.X = X
        self.top = top
        self.bottom = bottom

    def __repr__(self):
        X = self.X
        top = self.top
        bottom = self.bottom
        items = []
        if isinstance(X, dict):
            o = "{"
            c = "}"
            keys = sorted(list(X.keys()))
            for k in keys:
                items.append("%s: %s" % (k, X[k]))
        elif isinstance(X, list):
            o = "["
            c = "]"
            items = [str(x) for x in X]
        if len(items) > top + bottom:
            c += "\n(%s lines in total)" % len(items)
            items = items[:top] + ["..."] + items[-bottom:]
        return o + ",\n ".join(items) + c

