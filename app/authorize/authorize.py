from flask import Blueprint, redirect, render_template, flash, request, url_for,session,abort,jsonify
from flask_login import login_required, logout_user, current_user, login_user
from .. import login_manager
from .. import db,mail
from .forms import *
from .models import User
from flask import current_app as app
import uuid
from flask_mail import Message
from werkzeug.security import check_password_hash, generate_password_hash
import requests
from flask_jwt_extended import create_access_token


authorize_bp = Blueprint("authorize_bp", __name__, template_folder="templates", static_folder="static")

@authorize_bp.route('/')
def index():
    return redirect(url_for('authorize_bp.login'))

@authorize_bp.route("/signup", methods=["GET", "POST"])
def signup():
    """
    User sign-up page.

    GET requests serve sign-up page.
    POST requests validate form & user creation.
    """
    form = SignupForm()
    if form.validate_on_submit():
        # secret_response = request.form['g-recaptcha-response']
        # verify_response=requests.post(url=f'{app.config["VERIFY_URL"]}?secret={app.config["SECRET_KEY"]}&response={secret_response}').json()
        # print(verify_response)
        # if verify_response['success']==False or verify_response['score'] < 0.5:
        #     abort(401)
        # else:
        existing_user = User.query.filter_by(email=form.email.data).first()
        if existing_user is None:
            user = User(
                name=form.name.data,last_name=form.last_name.data, email=form.email.data,user_type="2"
            )
            user.set_password(form.password.data)
            verificacion_email = str(uuid.uuid4())
            user.set_verification_email(verificacion_email=verificacion_email)
            db.session.add(user)
            db.session.commit()  # Create new user
            #login_user(user)  # Log in as newly created user
            email_info = Message('Activate your account', sender = app.config['MAIL_USERNAME'], recipients = [form.email.data])
            activate_link = app.config['DOMAIN']+ url_for('authorize_bp.activate', email=str(form.email.data), code=verificacion_email)
            email_info.body = render_template('authorize/activation-email-template.html', link=activate_link)
            email_info.html = render_template('authorize/activation-email-template.html', link=activate_link)
            mail.send(email_info)
            return render_template("authorize/validar_email.html")
        flash("Este correo ya esta asignado a una cuenta")
    return render_template("authorize/signup.html", title="Create an Account.", form=form, template="signup-page", body="Sign up for a user account.",site_key=app.config["SITE_KEY"],)

@authorize_bp.route('/login/activate/<string:email>/<string:code>', methods=['GET'])
def activate(email, code):
    form = LoginForm()
    user = User.query.filter_by(email=email,verificacion_email=code).first()
    user.set_activated_account(activated_account="1")
    db.session.add(user)
    db.session.commit()
    #if user and user.chek_verificacion_email(code=code):
    session['loggedin'] = True
    session['nombre']=user.name
    session['id'] = user.id
    session['user_type'] = user.user_type
    login_user(user)
    flash("Tu cuenta ha sido activada, Bienvenido a Matchealo")
    return redirect(url_for("inicio_bp.inicio"))

@authorize_bp.route("/login", methods=["GET", "POST"])
def login():
    """
    Log-in page for registered users.
    GET requests serve Log-in page.
    POST requests validate and redirect user to dashboard.
    """
    # Bypass if user is logged in
    if current_user.is_authenticated:
        return redirect(url_for("inicio_bp.inicio"))
    form = LoginForm()
    # Validate login attempt
    if form.validate_on_submit():
        user = User.query.filter_by(email=form.email.data).first()
        print(user)
        if user and user.check_password(password=form.password.data):
            if user.activated_account=='1':
                login_user(user)
                next_page = request.args.get("next")
                session['loggedin'] = True
                session['nombre']=user.name
                session['id'] = user.id
                session['user_type'] = user.user_type
                session['admin'] = 'No'
                return redirect(url_for("inicio_bp.inicio"))
            if user.activated_account=='0':
                login_user(user)
                next_page = request.args.get("next")
                return redirect(url_for("authorize_bp.login"))
        flash("Invalid username/password combination")
        return redirect(url_for("authorize_bp.login"))
    return render_template("authorize/login.html", form=form, title="Log in.", template="login-page", body="Log in with your User account.",)

@authorize_bp.route("/api_login", methods=["POST"])
def login2():
    email = request.json.get("email", None)
    password = request.json.get("password", None)
    #user = User.query.filter_by(username=username).one_or_none()
    user = User.query.filter_by(email=email).one_or_none()
    if not user or not user.check_password(password):
        return jsonify("Wrong username or password"), 401
    login_user(user)
    # Notice that we are passing in the actual sqlalchemy user object here
    access_token = create_access_token(identity=current_user.email)
    return jsonify(access_token=access_token)

@authorize_bp.route('/logout')
@login_required
def logout():
    """User log-out logic."""
    logout_user()
    session['loggedin'] = False
    session['nombre']=""
    session['id'] = ""
    session['user_type'] = ""
    return redirect(url_for('authorize_bp.login'))

@login_manager.user_loader
def load_user(user_id):
    """Check if user is logged-in upon page load."""
    if user_id is not None:
        return User.query.get(user_id)
    return None

@login_manager.unauthorized_handler
def unauthorized():
    """Redirect unauthorized users to Login page."""
    flash("You must be logged in to view that page.")
    return redirect(url_for("authorize_bp.login"))

# -------proceso de recuperacion de contraseña
@authorize_bp.route("/recuperar", methods=["GET", "POST"])
def recuperar():
    form = CorreoForm()
    if form.validate_on_submit():
        existing_email = User.query.filter_by(email=form.email.data).first()
        if existing_email is not None:
            #login_user(user)  # Log in as newly created user
            email_info = Message('Cambiar contraseña', sender = app.config['MAIL_USERNAME'], recipients = [form.email.data])
            activate_link = app.config['DOMAIN']+ url_for('authorize_bp.password', email=str(form.email.data))
            email_info.body = render_template('authorize/activation-email-template.html', link=activate_link)
            email_info.html = render_template('authorize/activation-email-template.html', link=activate_link)
            mail.send(email_info)
            return render_template("authorize/validar_email.html")
            #return redirect(url_for("dashboard_bp.dashboard"))
        flash("correo no se encuentra registrado en el .")
    return render_template(
        "authorize/correo_recuperacion.html",
        title="Create an Account.",
        form=form,
        template="signup-page",
        body="Sign up for a user account.",
    )

@authorize_bp.route('/password/<string:email>', methods=['GET'])
def password(email):
    form = PassForm()
    user = User.query.filter_by(email=email).first()
    return render_template(
        "authorize/password.html",
        form=form,
        user=str(user.email),
        title="Create an Account.",
        template="signup-page",
        body="Sign up for a user account.",
    )

@authorize_bp.route('/word/<string:email>', methods=['POST'])
def word(email):
    form = PassForm()
    user = User.query.filter_by(email=email).first()
    if form.validate_on_submit():
        password=generate_password_hash(form.password.data, method="sha256")
        update=User.query.filter_by(email=user.email).update({User.password:password,})
        db.session.commit() 
        login_user(user)
        next_page = request.args.get("next")
        session['loggedin'] = True
        session['nombre']=user.name
        session['id'] = user.id
        session['user_type'] = user.user_type
        return redirect(next_page or url_for("explorar_bp.explorar",ix=1))
    return render_template(
        "authorize/password.html",
        title="Create an Account.",
        user=str(user.email),
        form=form,
        template="signup-page",
        body="Sign up for a user account.",
    )
