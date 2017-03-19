import configparser
import sqlite3
import os
import smtplib
import email
import uuid

config = configparser.ConfigParser()
config.read("/users/mathreview/settings.ini")

Email = config['Global']['Email']
Term = config['Global']['Term']
Database = sqlite3.connect("/users/mathreview/database/{0}.sqlite".format(Term))

# Initial database setup
_temp_cursor = Database.cursor()
_temp_cursor.execute('''CREATE TABLE IF NOT EXISTS submissions
	(id INTEGER PRIMARY KEY AUTOINCREMENT,
	 author TEXT,
	 email TEXT,
	 secret TEXT,
         title TEXT,
         abstract TEXT,
	 source BLOB,
	 pdf BLOB)''')
Database.commit()

def make_disk_copy(id, author, email, secret, title, abstract, source, pdf):
  target = "/users/mathreview/submissions/{0}/{1}".format(Term, id)
  os.mkdir(target)

  description = ''' Submission #{0}:
From: {1} <{2}>

{3}

Abstract
====================================
{4}
'''.format(id, author, email, title, abstract)
  
  with file(target + "/description.txt") as f:
    f.write(description)
  
  with file(target + "/source.tex") as f:
    f.write(source)

  with file(target + "/contents.pdf") as f:
    f.write(pdf)

def send_email(id, author, email, secret, title, abstract):
  contents = '''Dear {0},

Your submission to the Waterloo Math Review, {1}, has been successfully received
by our server.  Below find a copy of your abstract as submitted.  In case of
errors or omissions, please resubmit your paper.  We will review the latest
version of this paper that you submit.

Submission #{5}: {1}

Abstract
====================================
{2}
'''.format(author, title, abstract, id)
  message = email.MIMEMultipart()
  message['From'] = "Waterloo Math Review <{0}>".format(Email)
  message['To'] = "{0} <{1}>".format(author, email)
  message['Subject'] = "Submission #{0} - {1}".format(id, title)
  message.attach(email.MIMEText(contents, "plain"))

  server = smtplib.SMTP("localhost", 25)
  server.sendmail(Email, email, message.as_string())


# New submission function: (author, email, title, abstract, source, pdf) --> (id, secret)
def new_submission(author, email, title, abstract, source, pdf):
  cursor = Database.cursor()
  secret = uuid.uuid1().hex
  cursor.execute('''INSERT INTO submissions (author, email, secret, title, abstract, source, pdf) VALUES (?,?,?,?,?,?,?)''',
		 (author, email, secret, title, abstract, source, pdf))
  Database.commit()
  id = cursor.lastrowid
  
  # Save a copy
  make_disk_copy(id, author, email, secret, title, source, pdf)
  send_email(id, author, email, secret, title, source, pdf) 
