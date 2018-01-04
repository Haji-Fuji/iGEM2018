import flask
import pandas as pd
from sklearn.linear_model import LinearRegression

train_X = pd.read_csv("blood_fat.csv")
train_y = train_X.pop("blood fat")
model = LinearRegression().fit(train_X, train_y)


app = flask.Flask(__name__)


@app.route("/")
def index():
    return app.send_static_file("index.html")


@app.route("/bloodfat")
def bloodfat():
    age = flask.request.args.get("age", default=None, type=int)
    weight = flask.request.args.get("weight", default=None, type=float)
    if age is None or weight is None:
        return flask.jsonify({
            "code": 400,
            "msg": "Bad Request"
        })

    x = [age, weight]
    blood_fat = model.predict(x)[0]
    return flask.jsonify({
        "code": 200,
        "msg": "OK",
        "result": blood_fat
    })


if __name__ == "__main__":
    app.run(debug=True)
