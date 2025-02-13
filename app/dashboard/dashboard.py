"""Instantiate a Dash app."""
import dash
from dash import dash_table
from dash import dcc
from dash import html

import plotly.io as pio
import pandas as pd

from .data import create_dataframe

from flask import Blueprint, redirect, render_template, url_for,request,session,jsonify
from flask_login import current_user, login_required, logout_user, login_required
from .. import login_manager
from .. import db
from flask import current_app as app

from sklearn import cluster
import networkx as nx
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib import cm,colors
import mpld3
import seaborn as sns
import pandas as pd
import numpy as np
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score

import matplotlib
# https://matplotlib.org/stable/users/explain/backends.html
matplotlib.use('Agg')

dashboard_bp = Blueprint("dashboard_bp", __name__, template_folder="templates", static_folder="static")

@dashboard_bp.route("/dashboard", methods=["GET", "POST"])
def dashboard():
    # Obtener los datos y crear el dataframe df
    df = create_dataframe()
    # Renderizar el template con los datos necesarios
    plot_div = pio.to_html(
        {
            "data": [
                {
                    "x": df["tiempo"],
                    "y": df["distancia_promedio"],
                    "text": df["frecuencia_cardiaca"],
                    "mode": "markers",
                    "name": "Actividad Física",
                    "type": "scatter",
                }
            ],
            "layout": {
                "height": 700,
                "margin": {
                    "pad": 50
                },
                "xaxis": {"title": "Tiempo"},
                "yaxis": {"title": "Distancia Promedio"},
            },
        },
        full_html=False
    )
    return render_template('dashboard/dashboard.html', plot_div=plot_div, table=df.to_dict(orient='records')) 

@app.route('/data-source')
def data_source():
    df = create_dataframe()
    data = df.to_dict(orient='records')
    return jsonify(data)


def draw_communities(G, membership, pos):
    fig, ax = plt.subplots(figsize=(7, 3))
    club_dict = defaultdict(list)
    for student, club in enumerate(membership):
        club_dict[club].append(student)
    norm = colors.Normalize(vmin=0, vmax=len(club_dict.keys()))
    for club, members in club_dict.items():
        nx.draw_networkx_nodes(G, pos,
                            nodelist=members,
                            node_color=[cm.jet(norm(club))],
                            node_size=50,
                            alpha=0.8,
                            ax=ax)   
    plt.title("Zachary's Karate Club")
    nx.draw_networkx_edges(G, pos, alpha=0.5, ax=ax)


@dashboard_bp.route('/karate_graph')
def karate_graph():
    G = nx.karate_club_graph()
    pos = nx.spring_layout(G)
    y_true = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    draw_communities(G, y_true, pos)
    plt.tight_layout()
    graph_html = mpld3.fig_to_html(plt.gcf())
    return render_template('dashboard/karate.html', graph_html=graph_html)

