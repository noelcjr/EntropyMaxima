"""Sign-up & log-in forms."""
from flask_wtf import FlaskForm
from wtforms import PasswordField, StringField, SubmitField
from wtforms.validators import DataRequired, Email, EqualTo, Length, Optional

class SignupForm(FlaskForm):
    """User Sign-up Form."""
    name = StringField("First name", validators=[DataRequired()])
    last_name = StringField("Last name", validators=[DataRequired()])
    email = StringField(
        "Email",
        validators=[
            Length(min=6),
            Email(message="Enter a valid email."),
            DataRequired(),
        ],
    )
    password = PasswordField(
        "Password",
        validators=[
            DataRequired(),
            Length(min=6, message="Select a stronger password."),
        ],
    )
    confirm = PasswordField(
        "Confirm password",
        validators=[
            DataRequired(),
            EqualTo("password", message="Passwords must match."),
        ],
    )
    submit = SubmitField("Register")

class LoginForm(FlaskForm):
    """User Log-in Form."""
    
    email = StringField(
        "Email", validators=[DataRequired(), Email(message="Enter a valid email.")]
    )
    password = PasswordField("Password", validators=[DataRequired()])
     
    submit = SubmitField("Log in")



class PassForm(FlaskForm):
    """User Log-in Form."""
 
    password = PasswordField(
        "Password",
        validators=[
            DataRequired(),
            Length(min=6, message="Enter a secure password"),
        ],
    )
    confirm = PasswordField(
    "Confirm password",
    validators=[
        DataRequired(),
        EqualTo("password", message="Passwords do not match"),
    ],
)


    submit = SubmitField("Cambiar")

class CorreoForm(FlaskForm):
    """User Log-in Form."""
 
    email = StringField(
        "Email", validators=[DataRequired(), Email(message="Enter a valid email.")]
    )     
    submit = SubmitField("Enviar")
