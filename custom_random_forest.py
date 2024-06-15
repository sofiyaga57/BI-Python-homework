from concurrent.futures import ProcessPoolExecutor
import numpy as np
from sklearn.tree import DecisionTreeClassifier
from sklearn.base import BaseEstimator
from functools import partial


class RandomForestClassifierCustom(BaseEstimator):
    def __init__(self, n_estimators=10, max_depth=None, max_features=None, random_state=144):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state
        self.trees = []
        self.feat_ids_by_tree = []

    def _fit_loop(self, i, X, y):
        np.random.seed(self.random_state + i)
        feat_ids = np.random.choice(range(X.shape[1]), self.max_features, replace=False)
        x_ids_bootstrap = np.random.choice(range(X.shape[0]), X.shape[0], replace=True)

        clf = DecisionTreeClassifier(max_depth=self.max_depth, max_features=self.max_features,
                                     random_state=self.random_state)
        clf.fit(X[x_ids_bootstrap, :][:, feat_ids], y[x_ids_bootstrap])
        return clf, feat_ids

    def fit(self, X, y, n_jobs=1):
        self.classes_ = sorted(np.unique(y))
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            results = list(executor.map(partial(self._fit_loop, X=X, y=y), range(self.n_estimators)))
        self.trees, self.feat_ids_by_tree = zip(*results)
        return self

    def _predict_proba_loop(self, X, clf, feat_ids):
        return clf.predict_proba(X[:, feat_ids])

    def predict_proba(self, X, n_jobs=1):
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            futures = [executor.submit(self._predict_proba_loop, X, self.trees[i], self.feat_ids_by_tree[i]) for i in
                       range(self.n_estimators)]
            y_pred_s = [future.result() for future in futures]
        return np.mean(y_pred_s, axis=0)

    def predict(self, X, n_jobs=1):
        probas = self.predict_proba(X, n_jobs)
        predictions = np.argmax(probas, axis=1)
        return predictions
