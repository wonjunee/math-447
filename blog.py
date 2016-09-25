import os
import re
import random
import hashlib
import hmac
from string import letters

import webapp2
import jinja2

from google.appengine.ext import db

template_dir = os.path.join(os.path.dirname(__file__), 'templates')
jinja_env = jinja2.Environment(loader = jinja2.FileSystemLoader(template_dir),
                               autoescape = True)

secret = 'helloworld'

def render_str(template, **params):
    t = jinja_env.get_template(template)
    return t.render(params)

class BlogHandler(webapp2.RequestHandler):
    def write(self, *a, **kw):
        self.response.out.write(*a, **kw)

    def render_str(self, template, **params):
        return render_str(template, **params)

    def render(self, template, **kw):
        self.write(self.render_str(template, **kw))

# Home
class Home(BlogHandler):
    def get(self):
        self.render('home.html')

# Project Pages
class ProjectPage(BlogHandler):
    def get(self,project_num):
        self.render('project{}.html'.format(project_num))

# Codes Pages
class CodesPage(BlogHandler):
    def get(self, project_num, code_page):
        self.render('{}.html'.format(code_page))

app = webapp2.WSGIApplication([
                               ('/', Home),
                               ('/project/([0-9])', ProjectPage),
                               ('/project/([0-9])/([0-9a-zA-z]+)', CodesPage),
                               ],
                              debug=True)
