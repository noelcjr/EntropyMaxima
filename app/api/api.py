from flask import Flask, Blueprint, request, jsonify
from flask import send_from_directory
import os
from .models import Calls
from flask_login import current_user, login_required, logout_user, login_required
from .. import login_manager
from ..authorize.models import User
from .. import db
from flask_jwt_extended import JWTManager, create_access_token, jwt_required, get_jwt_identity

api_bp = Blueprint('api_bp', __name__, template_folder='templates', static_folder='static')

def validate_call(id):
    try:
        item = Calls.query.get(id)
        del item.__dict__['_sa_instance_state']
        return True
    except:
        return False

@api_bp.route('/static/todos.json')
def send_json():
    return send_from_directory('static', 'todos.json')

@api_bp.route('/get/<id>', methods=['GET'])
@jwt_required()
def get_item(id):
    current_user_email = get_jwt_identity()
    user = User.query.filter_by(email=current_user_email).first()
    try:
        item = Calls.query.filter_by(id=id, user_id=user.id).first()
        del item.__dict__['_sa_instance_state']
        return jsonify(item.__dict__)
    except:
        response = jsonify({'error': 'Resource not found'})
        response.status_code = 404
        return response

# Lo quite para evitar que nos bajen toda la base de datos y nos roben informacion.
# por las mismas razones debemos no permitir mas de un numero de GETs por dia.
@api_bp.route('/get_all', methods=['GET'])
@jwt_required()
def get_items():
    items = []
    current_user_email = get_jwt_identity()
    user = User.query.filter_by(email=current_user_email).first()
    try:
        for item in db.session.query(Calls).filter(Calls.user_id == user.id).all():
            del item.__dict__['_sa_instance_state']
            items.append(item.__dict__)
        return jsonify(items)
    except:
        response = jsonify({'error': 'Resource not found'})
        response.status_code = 404
        return response

@api_bp.route('/post', methods=['POST'])
@jwt_required()
def create_item():
    body = request.get_json()
    current_user_email = get_jwt_identity()
    user = User.query.filter_by(email=current_user_email).first()
    try:
        db.session.add(Calls(key=body['key'],user_id = user.id, created=body['created'],agency= body['agency'],agency_name= body['agency_name'],complaint_type= body['complaint_type'], descriptor=body['descriptor'], incident_zip=body['incident_zip'],city= body['city'], latitud=body['latitud'],longitude= body['longitude']))
        db.session.commit()
        response = jsonify({'success': 'Item created'})
        response.status_code = 200
        return response
    except:
        response = jsonify({'error': 'Conflict, key already exists'})
        response.status_code = 409
        return response

@api_bp.route('/put/<id>', methods=['PUT'])
@jwt_required()
def update_item(id):
    body = request.get_json()
    current_user_email = get_jwt_identity()
    user = User.query.filter_by(email=current_user_email).first()
    if validate_call(id)==True:
        try:
            db.session.query(Calls).filter_by(id=id, user_id= user.id).update(
                dict(key=body['key'], created=body['created'],agency= body['agency'],agency_name= body['agency_name'],complaint_type= body['complaint_type'], descriptor=body['descriptor'], incident_zip=body['incident_zip'],city= body['city'], latitud=body['latitud'],longitude= body['longitude']))
            db.session.commit()
            response = jsonify({'success': 'item updated'})
            response.status_code = 200
            return response
        except:
            return "Error"
    else: 
        response = jsonify({'error': 'Resource not found'})
        response.status_code = 404
        return response
    
@api_bp.route('/delete/<id>', methods=['DELETE'])
@jwt_required()
def delete_item(id):
    current_user_email = get_jwt_identity()
    user = User.query.filter_by(email=current_user_email).first()
    if validate_call(id)==True:
        try:
            db.session.query(Calls).filter_by(id=id,user_id=user.id).delete()
            db.session.commit()
            response = jsonify({'success': 'item deleted'})
            response.status_code = 204
            return response
        except:
            return "Error"
    else:
        response = jsonify({'error': 'Resource not found'})
        response.status_code = 404
        return response


