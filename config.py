"""Flask config."""
from os import environ, path
import os
import getpass
from dotenv import load_dotenv

BASE_DIR = path.abspath(path.dirname(__file__))
load_dotenv(path.join(BASE_DIR, ".env"))

class Config:
    """Flask configuration variables."""

    # General Config
    FLASK_APP = "wsgi.py"
    FLASK_ENV = environ.get("FLASK_ENV")
    SECRET_KEY = environ.get("SECRET_KEY")

    # Assets
    LESS_BIN = environ.get("LESS_BIN")
    ASSETS_DEBUG = environ.get("ASSETS_DEBUG")
    LESS_RUN_IN_DEBUG = environ.get("LESS_RUN_IN_DEBUG")

    # Static Assets
    STATIC_FOLDER = "static"
    TEMPLATES_FOLDER = "templates"
    COMPRESSOR_DEBUG = environ.get("COMPRESSOR_DEBUG")

    usr = getpass.getuser()
    PROJECT_ID = 'Entropy Maxima'
    if usr == 'noel': 
        SECRET_KEY = os.environ.get('secrete_key')
        SITE_KEY = os.environ["site_key"]
        MYSQL_HOST = '127.0.0.1'
        MYSQL_USER = 'noel'
        MYSQL_PASSWORD = os.environ.get('cloudSQL_password_em')
        DOMAIN ='http://localhost:5000'
        SQLALCHEMY_DATABASE_URI = "mysql+pymysql://noel:"+MYSQL_PASSWORD+"@127.0.0.1:5434/EntropyMaxima"
        SQLALCHEMY_TRACK_MODIFICATIONS = False
        MAIL_SERVER = 'smtp.gmail.com'
        MAIL_PORT = 465
        MAIL_USERNAME = 'contact@xmatchplus.com'
        MAIL_PASSWORD = os.environ['MAIL_PASSWORD_XMATCHPLUS']
        MAIL_USE_TLS = False
        MAIL_USE_SSL = True
        GITHUB_ID=os.environ.get('GITHUB_ID')
        GITHUB_SECRET=os.environ.get('GITHUB_SECRET')
        OAUTHLIB_INSECURE_TRANSPORT=1
    elif usr == 'libardo':
        SECRET_KEY = os.environ.get('secrete_key')
        SITE_KEY = os.environ["site_key"]
        MYSQL_HOST = '127.0.0.1' #'db'
        MYSQL_USER = 'permutas'  #'root'
        MYSQL_PASSWORD = os.environ.get('cloudSQL_password2')
        DOMAIN ='http://localhost:5000'
        SQLALCHEMY_DATABASE_URI = "mysql://libardo:"+MYSQL_PASSWORD+"@127.0.0.1:5437/ritmotion"
        SQLALCHEMY_TRACK_MODIFICATIONS = False
        MAIL_SERVER = 'smtp.gmail.com'
        MAIL_PORT = 465
        MAIL_USERNAME = 'admin@matchealo.com'
        MAIL_PASSWORD = os.environ['mail_password']
        MAIL_USE_TLS = False
        MAIL_USE_SSL = True
    elif usr == 'root':
        SECRET_KEY = os.environ["SECRET_KEY"]
        SITE_KEY = os.environ["SITE_KEY"]
        VERIFY_URL= 'https://www.google.com/recaptcha/api/siteverify'
        MYSQL_HOST = os.environ["PERMUTAS_DB_HOST"]
        MYSQL_USER = os.environ["PERMUTAS_DB_USER"]
        MYSQL_PASSWORD = os.environ["PERMUTAS_DB_PASSWORD"]
        BUCKET_NAME = os.environ["PERMUTAS_BUCKET"]
        DOMAIN ='https://ritmomotion.com'
        UPLOAD_FOLDER = '/app/matchealo/static/fotos/'
        SQLALCHEMY_DATABASE_URI = "mysql+pymysql://noel:"+MYSQL_PASSWORD+"@127.0.0.1:3306/ritmotion"
        SQLALCHEMY_TRACK_MODIFICATION = False
        ALLOWED_EXTENSIONS = ['.png', '.jpg','.jpeg']
        MAX_CONTENT_LENGTH = 4096 * 4096
        MAIL_SERVER = 'smtp.gmail.com'
        MAIL_PORT = 465
        MAIL_USERNAME = 'admin@matchealo.com'
        MAIL_PASSWORD = os.environ['ADMIN_GMAIL_PASSWORD']
        MAIL_USE_TLS = False
        MAIL_USE_SSL = True
    elif usr == 'jhon':
        SECRET_KEY = os.environ.get('agenda_crm_key')
        SITE_KEY = os.environ["site_key"]
        VERIFY_URL= 'https://www.google.com/recaptcha/api/siteverify'
        MYSQL_HOST = '127.0.0.1' #'db'
        MYSQL_USER = 'permutas'  #'root'
        MYSQL_PASSWORD = os.environ.get('cloudSQL_password2')
        DOMAIN ='http://localhost:5000'
        SQLALCHEMY_DATABASE_URI = "mysql://jhon:"+MYSQL_PASSWORD+"@127.0.0.1:5434/ritmotion"
        SQLALCHEMY_TRACK_MODIFICATIONS = False
        MAIL_SERVER = 'smtp.gmail.com'
        MAIL_PORT = 465
        MAIL_USERNAME = 'admin@matchealo.com'
        MAIL_PASSWORD = os.environ['mail_password']
        MAIL_USE_TLS = False
        MAIL_USE_SSL = True
    else:
        print('ERROR: No valid value for PROJECT_ENV')
    SQLALCHEMY_TRACK_MODIFICATIONS = False 
