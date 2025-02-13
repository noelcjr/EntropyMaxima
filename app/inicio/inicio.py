from flask import Blueprint, redirect, render_template, url_for,request,session,jsonify
from flask_login import current_user, login_required, logout_user, login_required
from .. import login_manager
from .forms import PerfilForm
# from ..authorize import models
#from .models import Persona
from .. import db
from flask import current_app as app
import os
import base64
from ..authorize.models import User
from flask_jwt_extended import JWTManager, create_access_token, jwt_required, get_jwt_identity

inicio_bp = Blueprint("inicio_bp", __name__, template_folder="templates", static_folder="static")

@inicio_bp.route("/inicio", methods=["GET", "POST"])
@login_required
def inicio():
    #  return render_template(
    #  "inicio.html",
    #  current_user=current_user,
    #  body="inicio.",
    # #  foto_perfil=foto_perfil
    #  )
    if session["admin"]=="No":
        return redirect(url_for("dashboard_bp.dashboard"))
    elif session["admin"]=="Si":
        return render_template(
        "page_blank.html",
        )

@inicio_bp.route("/api_key", methods=["GET", "POST"])
@login_required
def api_key():
     return render_template(
     "obtener_api.html",
        api_key=current_user.token
     )

@inicio_bp.route("/update_api_key", methods=["POST"])
@login_required
def update_api_key():
    access_token = create_access_token(identity=current_user.email)
    update = User.query.filter_by(email=current_user.email).update({User.token: access_token})
    db.session.commit()
    return jsonify({"api_key": access_token}), 200

@inicio_bp.route("/inicio_admin", methods=["GET", "POST"])
@login_required
def inicio_admin():
    """ inicio """
    session['admin']="Si"
    return redirect(url_for("inicio_bp.inicio"))

@inicio_bp.route("/quitar_admin", methods=["GET", "POST"])
@login_required
def quitar_admin():
    """ inicio """
    session['admin']="No"
    return redirect(url_for("inicio_bp.inicio"))

