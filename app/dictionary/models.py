from .. import db
from flask_login import UserMixin
from sqlalchemy.orm import relationship
from sqlalchemy import ForeignKey

# _ = db.Column(db., nullable=False, unique=True)
# _id = db.Column(ForeignKey(""), primary_key=True)
# _ = relationship("", back_populates="")

class Dict_ST_Association(UserMixin, db.Model):
    """ Associates Dictionarys to Semantic Types """
    __tablename__ = "dict_st_association"
    dict_id = db.Column(ForeignKey("dictionary.id"), primary_key=True)
    st_id = db.Column(ForeignKey("semantic_type.id"), primary_key=True)
    date_created = db.Column(db.String(12))
    dictionary = relationship("Dictionary", back_populates="semantic_types")
    semantic_type = relationship("Semantic_Type", back_populates="dicitonaries")

class Dictionary(UserMixin, db.Model):
    """ Fro reference adn custom dictionary, metadata only"""
    __tablename__ = "dicitonary"
    id = db.Column(db.Integer, nullable=False, primary_key=True)
    dicitonary_Name = db.Column(db.String(25), nullable=False, unique=False)
    dictionary_trunk = db.Column(db.String(25), nullable=False, unique=False)
    dictioanry_branch = db.Column(db.String(25), nullable=False, unique=False)
    owner = db.Column(db.String(50), nullable=False, unique=False)
    email = db.Column(db.String(50), nullable=False, unique=False)
    created_on = db.Column(db.String(12), index=False, nullable=False, unique=False)
    last_modified = db.Column(db.String(12), index=False, nullable=False, unique=True)
    created_on = db.Column(db.String(12), index=False, nullable=False, unique=False)
    description = db.Column(db.Text, nullable=True, unique=False)
    reference = db.Column(db.Integer, nullable=False, unique=False)
    semantic_types = relationship("Dict_ST_Association", back_populates="dictioanary")
    users = relationship("User_Dictionary_Association", back_populates="dictionary")
    terms_newterms = relationship("Term_NewTerm_Association", back_populates="dictionary")

class ST_Term_Association(UserMixin, db.Model):
    """ Associates Semantic Types and Terms from separate table"""
    __tablename__ = "st_term_association"
    st_id = db.Column(ForeignKey("semantic_type.id"), primary_key=True)
    term_id = db.Column(ForeignKey("term.id"), primary_key=True)
    date_created = db.Column(db.String(12), index=False, nullable=False, unique=False)
    semantic_type = relationship("Semantic_Type", back_populates="terms")
    term = relationship("Term", back_populates="semantic_types")

class Semantic_Type(UserMixin, db.Model):
    __tablename__ = "semantic_type"
    dict_id = db.Column(db.Integer, primary_key=True)
    SemType_Name = db.Column(db.String(100), nullable=False, unique=True)
    created_on = db.Column(db.String(12), index=False, nullable=False, unique=False)
    last_modified = db.Column(db.String(12), index=False, nullable=False, unique=False)
    description = db.Column(db.Text, nullable=True, unique=False)
    dictionaries = relationship("Dict_ST_Association", back_populates="semantic_type")
    terms = relationship("ST_Term_Association", back_populates="semantic_type")
    terms_newterms = relationship("Term_NewTerm_Association", back_populates="semantic_type")

class Term(UserMixin, db.Model):
    __tablename__ = "term"
    id = db.Column(db.Integer, primary_key=True)
    term = db.Column(db.String(200), nullable=False)
    code = db.Column(db.String(20), nullable=True)
    parent_code = db.Column(db.String(20), nullable=True)
    definition = db.Column(db.Text, nullable=True)
    display_name = db.Column(db.String(200), nullable=True)
    concept_status = db.Column(dbString(200)., nullable=True)
    concept_in_subset = db.Column(db.String(200), nullable=True)
    semantic_types = relationship("ST_Term_Association", back_populates="term")
    terms = relationship("Term_NewTerm_Association", back_populates="term")

class Synonym(UserMixin, db.Model):
    __tablename__ = "synonyms"
    id = db.Column(db.Integer, primary_key=True)
    term_id = db.Column(ForeignKey("term.id"))
    synonyms = db.Column(db.String(200))

