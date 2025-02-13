import os
from flask import Blueprint, Flask, render_template, request, redirect, url_for
#from flask_sqlalchemy import SQLAlchemy
from .models import Tag, Post, Comment
from .. import db

# https://www.digitalocean.com/community/tutorials/how-to-use-one-to-many-database-relationships-with-flask-sqlalchemy
# https://www.digitalocean.com/community/tutorials/how-to-use-many-to-many-database-relationships-with-flask-sqlalchemy

blog_bp = Blueprint("blog_bp", __name__, template_folder="templates", static_folder="static")

@blog_bp.route('/blog')
def index():
    posts = Post.query.all()
    return render_template('blog/index.html', posts=posts)

@blog_bp.route('/<int:post_id>/', methods=('GET', 'POST'))
def post(post_id):
    post = Post.query.get_or_404(post_id)
    if request.method == 'POST':
        comment = Comment(content=request.form['content'], post=post)
        db.session.add(comment)
        db.session.commit()
        return redirect(url_for('blog_bp.post', post_id=post.id))
    return render_template('blog/post.html', post=post)

@blog_bp.route('/comments/')
def comments():
    comments = Comment.query.order_by(Comment.id.desc()).all()
    return render_template('blog/comments.html', comments=comments)

@blog_bp.route('/comments/<int:comment_id>/delete')
def delete_comment(comment_id):
    comment = Comment.query.get_or_404(comment_id)
    post_id = comment.post.id
    db.session.delete(comment)
    db.session.commit()
    return redirect(url_for('blog_bp.post', post_id=post_id))

@blog_bp.route('/tags/<tag_name>/')
def tag(tag_name):
    tag = Tag.query.filter_by(name=tag_name).first_or_404()
    return render_template('blog/tag.html', tag=tag)

@blog_bp.route('/load_db', methods=('GET', 'POST'))
def load_db():
    #db.drop_all()
    #db.create_all()

    post1 = Post(title='Post The First', content='Content for the first post')
    post2 = Post(title='Post The Second', content='Content for the Second post')
    post3 = Post(title='Post The Third', content='Content for the third post')
    life_death_post = Post(title='A post on life and death', content='life and death')
    joy_post = Post(title='A post on joy', content='joy')

    comment1 = Comment(content='Comment for the first post', post=post1)
    comment2 = Comment(content='Comment for the second post', post=post2)
    comment3 = Comment(content='Another comment for the second post', post_id=2)
    comment4 = Comment(content='Another comment for the first post', post_id=1)

    tag1 = Tag(name='animals')
    tag2 = Tag(name='tech')
    tag3 = Tag(name='cooking')
    tag4 = Tag(name='writing')
    life_tag = Tag(name='life')
    death_tag = Tag(name='death')
    joy_tag = Tag(name='joy')
    
    post1.tags.append(tag1)  # Tag the first post with 'animals'
    post1.tags.append(tag4)  # Tag the first post with 'writing'
    post3.tags.append(tag3)  # Tag the third post with 'cooking'
    post3.tags.append(tag2)  # Tag the third post with 'tech'
    post3.tags.append(tag4)  # Tag the third post with 'writing'
    life_death_post.tags.append(life_tag)
    life_death_post.tags.append(death_tag)
    joy_post.tags.append(joy_tag)

    db.session.add_all([post1, post2, post3, life_death_post, joy_post])
    db.session.add_all([comment1, comment2, comment3, comment4])
    db.session.add_all([tag1, tag2, tag3, tag4, life_tag, death_tag, joy_tag])

    db.session.commit()
    posts = Post.query.all()
    return render_template('blog/index.html', posts=posts)
