from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, FloatField, validators, IntegerField, DateField, BooleanField, SelectField, DateTimeField
from wtforms.validators import DataRequired, Email, EqualTo, Length, Optional
# from wtforms.fields.html5 import DateField

from wtforms.widgets import TextArea

class PerfilForm(FlaskForm):
    nombre_1 = StringField("Primer Nombre:", validators=[DataRequired(), validators.Length(min=1, max=50)])
    nombre_2 = StringField("Segundo Nombre:", validators=[Optional(), validators.Length(min=1, max=50)])
    apellido_1 = StringField("Primer Apellido", validators=[DataRequired(), validators.Length(min=1, max=50)]) 
    apellido_2= StringField("Segundo Apellido", validators=[Optional(), validators.Length(min=1, max=50)]) 
    identificacion = StringField("Identificacion", validators=[DataRequired(), validators.Length(min=8, max=10)]) 
    estado_civil=SelectField(u'Estado Civil', choices=[('', 'Estado'),('soltero', 'soltero'), ('casado', 'casado')],validators=[DataRequired()])
    direccion=StringField("Direccion:")
    ubicacion=StringField("Barrio/Localidad/Zona:", validators=[DataRequired(), validators.Length(min=1, max=100)])
    movil_1=StringField("Movil 1:", validators=[DataRequired(),validators.Length(min=8, max=10)]) 
    movil_2=StringField("Movil 2:", validators=[Optional(),validators.Length(min=8, max=10)])
    email=StringField("Email:", validators=[DataRequired(), validators.Length(min=1, max=50),Email("This field requires a valid email address")])
    email_2=StringField("Email 2:", validators=[Optional(), validators.Length(min=1, max=50),Email("This field requires a valid email address")])

    submit = SubmitField("Actualizar")
