"""Initialize Flask app."""
from flask import Flask
from flask_assets import Environment
from flask_sqlalchemy import SQLAlchemy
from flask_login import LoginManager
import os
import getpass
from flask_mail import Mail
from flask_jwt_extended import JWTManager
#from flask_dance.contrib.github import github
#from flask_oauth import OAuth

db = SQLAlchemy()
login_manager = LoginManager()
mail = Mail()
#oauth = OAuth()

usr = getpass.getuser()
if usr == 'noel':
    os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = '/home/noel/maximusentropy1-31a8b44dc130.json'
if usr == 'libardo':
    os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = '/mnt/c/Users/bejar/OneDrive/Escritorio/trumelabs/trumelabs/maximusentropy1-098c6eb28f61.json'
if usr == 'jhon':
    os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = '/mnt/c/Users/jhon/Documents/keysNube/maximusentropy1-e671bd64b748.json'

def init_app():
    """Construct core Flask application with embedded Dash app."""
    app = Flask(__name__, instance_relative_config=False)
    app.config.from_object("config.Config")
    db.init_app(app)
    login_manager.init_app(app)
    mail = Mail(app)
    assets = Environment()
    assets.init_app(app)
    app.config['JWT_ACCESS_TOKEN_EXPIRES'] = 3600 # Tiempo de expiración en segundos (1 hora)
    jwt = JWTManager(app)

    with app.app_context():
        # Import parts of our core Flask app
        from . import routes
        from .assets import compile_static_assets
        from .authorize import authorize
        from .authorize import githubOAuth 
        from .inicio import inicio
        from .dashboard import dashboard
        from .api import api
        from .dashboard import init_dashboard
        from .blog import routes

        # Import Dash application
        app.register_blueprint(authorize.authorize_bp)
        #app.register_blueprint(githubOAuth.github_bp, url_prefix="/login")
        app.register_blueprint(githubOAuth.github_bp)
        app.register_blueprint(inicio.inicio_bp)
        #app = init_dashboard(app)
        app.register_blueprint(dashboard.dashboard_bp)
        app.register_blueprint(api.api_bp)
        app.register_blueprint(routes.blog_bp)
        # Compile static assets
        compile_static_assets(assets)
        db.create_all()
        #print(f"---{os.getenv('GITHUB_ID')}--")        
        #@app.route("/github")
        #def login2():
        #    if not github.authorized:
        #        return redirect(url_for("github.login"))
        #    res = github.get("/user")
        #    return f"You are @{res.json()['login']} on GitHub"
        return app

#@app.route("/ping")
#def ping():
#    return jsonify(ping="pong")
