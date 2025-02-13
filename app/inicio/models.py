from .. import db
from flask_login import UserMixin
from werkzeug.security import check_password_hash, generate_password_hash

class Persona(UserMixin, db.Model):
    """Persona """
    __tablename__ = "persona"
    id = db.Column(db.Integer, primary_key=True)
    nombre_1 = db.Column(db.String(50), nullable=False)
    nombre_2 = db.Column(db.String(50))
    apellido_1 = db.Column(db.String(50),  nullable=False)
    apellido_2 = db.Column(db.String(50))
    identificacion = db.Column(db.String(10), unique=True, nullable=False)
    pais_expedicion = db.Column(db.String(40), nullable=False)
    fecha_expedicion= db.Column(db.DateTime, nullable=False, unique=False)
    fecha_nacimiento = db.Column(db.DateTime, unique=False, nullable=False)
    direccion = db.Column(db.String(100))
    genero = db.Column(db.String(1), nullable=False)
    estado_civil = db.Column(db.String(1), nullable=False)
    ubicacion = db.Column(db.String(50))
    estrato = db.Column(db.String(1), nullable=False)
    movil_1 = db.Column(db.String(10), unique=False, nullable=False)
    movil_2 = db.Column(db.String(10), unique=False)
    email_1= db.Column(db.String(60), unique=True, nullable=False)
    email_2 = db.Column(db.String(60), unique=False)

    def __repr__(self):
        return "<cliente {}>".format(self.id)

