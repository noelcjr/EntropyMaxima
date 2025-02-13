"""Routes for parent Flask app."""
from flask import current_app as app
from flask import render_template
from .authorize import authorize

@app.route("/")
def login():
    return authorize.login()

