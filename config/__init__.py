# # config/__init__.py

# import os

# # Load environment variables
# from dotenv import load_dotenv

# # Load environment variables from .env file if it exists
# dotenv_path = os.path.join(os.path.dirname(__file__), '.env')
# if os.path.exists(dotenv_path):
#     load_dotenv(dotenv_path)

# # Set default configuration settings
# DATABASE_URI = os.getenv('DATABASE_URI', 'sqlite:///myapp.db')
# SECRET_KEY = os.getenv('SECRET_KEY', 'mysecretkey')
# DEBUG = os.getenv('DEBUG', False)
