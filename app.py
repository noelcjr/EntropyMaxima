"""Application entry point."""
# dsrfvd
from app import init_app
import getpass

usr = getpass.getuser()
app = init_app()

if __name__ == "__main__":
    app.run(host="0.0.0.0",debug=True)
