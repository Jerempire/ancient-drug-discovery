"""Send V2 campaign report via Gmail SMTP."""
import sys
try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

import smtplib
import os
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from dotenv import load_dotenv
from pathlib import Path

load_dotenv(Path('C:/Users/jmj2z/Projects/intelligence/substack-digest/.env'))
pw = os.environ.get('GMAIL_APP_PASSWORD', '')
if not pw:
    print('ERROR: No GMAIL_APP_PASSWORD')
    sys.exit(1)

recipient = 'jerempire99@gmail.com'

msg = MIMEMultipart('alternative')
msg['From'] = recipient
msg['To'] = recipient
msg['Subject'] = 'Ancient Drug Discovery: V2 Full Campaign Report - mpnn06 Lead (1.8x IRAP Selectivity)'

html = open(Path(__file__).parent / 'campaign_report.html').read()
text = 'V2 ERAP2 Inhibitor Campaign - mpnn06 lead: ERAP2=0.815, ERAP1=0.326 (2.5x), IRAP=0.454 (1.8x). See HTML for full tables.'

msg.attach(MIMEText(text, 'plain', 'utf-8'))
msg.attach(MIMEText(html, 'html', 'utf-8'))

with smtplib.SMTP('smtp.gmail.com', 587) as server:
    server.starttls()
    server.login(recipient, pw)
    server.sendmail(recipient, recipient, msg.as_string())
    print('Email sent successfully')
