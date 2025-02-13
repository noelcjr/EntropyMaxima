import os
from flask import Flask, Blueprint, jsonify, redirect, url_for
from flask_dance.contrib.github import make_github_blueprint
from flask_dance.contrib.github import github

#github_bp = make_github_blueprint(client_id=os.environ.get("GITHUB_ID"),
#                                  client_secret=os.environ.get("GITHUB_SECRET"),)

github_bp = Blueprint("github_bp", __name__)

@github_bp.route("/ping")
def ping():
    return jsonify(ping="pong")

#@github_bp.route("/github")
#def github_login():
#    if not github.authorized:
#        return redirect(url_for("github.login"))
#    res = github.get("/user")
#    return f"You are @{res.json()['login']} on GitHub"
