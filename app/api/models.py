from .. import db
from flask_login import UserMixin

class Calls(UserMixin, db.Model):
    id = db.Column(db.Integer, primary_key=True)
    key = db.Column(db.String(12), unique=True, nullable=False)
    user_id = db.Column(db.Integer, nullable=False)  # Columna para almacenar el ID del usuario
    created = db.Column(db.DateTime, nullable=False)
    agency = db.Column(db.String(10), unique=False, nullable=False)
    agency_name = db.Column(db.String(80), unique=False, nullable=False)
    complaint_type = db.Column(db.String(80), unique=False, nullable=False)
    descriptor = db.Column(db.String(200), unique=False, nullable=True)
    incident_zip =  db.Column(db.String(10), unique=False, nullable=False)
    city = db.Column(db.String(100), index=False, unique=False, nullable=False)
    latitud = db.Column(db.Float, nullable=False)
    longitude = db.Column(db.Float, nullable=False)

    def __repr__(self):
        return "<Calls {}>".format(self.id)