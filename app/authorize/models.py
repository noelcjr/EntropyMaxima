from .. import db
from flask_login import UserMixin
from werkzeug.security import check_password_hash, generate_password_hash

class User(UserMixin, db.Model):
    """User account model."""

    __tablename__ = "flasklogin-users"
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(100), nullable=False, unique=False)
    last_name = db.Column(db.String(40), )
    email = db.Column(db.String(40), unique=True, nullable=False)
    user_type= db.Column(db.String(40), unique=False, nullable=False)
    password = db.Column(db.String(200), primary_key=False, unique=False, nullable=False)
    verificacion_email = db.Column(db.String(255), unique=True)
    activated_account = db.Column(db.String(1), unique=False, nullable=False, default=False)
    token = db.Column(db.String(400),default=False)
    created_on = db.Column(db.DateTime, index=False, unique=False, nullable=True)
    last_login = db.Column(db.DateTime, index=False, unique=False, nullable=True)

    def set_password(self, password):
        """Create hashed password."""
        self.password = generate_password_hash(password, method="sha256")
    def check_password(self, password):
        """Check hashed password."""
        return check_password_hash(self.password, password)
    def check(self):
        """Check hashed password."""
        return self.password
    def set_verification_email(self, verificacion_email):
        """Crea una clave para verificar correos de usuario """
        self.verificacion_email = verificacion_email
    def set_activated_account(self, activated_account):
        self.activated_account = activated_account
    def chek_verificacion_email(self, code):
        if code == str(self.verificacion_email):
            return True
        else:
            return False
    def __repr__(self):
        return "<User {}>".format(self.name)
