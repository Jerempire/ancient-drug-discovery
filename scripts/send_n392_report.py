"""Send N392 selectivity screen HTML report via Gmail SMTP."""
import sys, os, smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from pathlib import Path
from dotenv import load_dotenv

try:
    sys.stdout.reconfigure(encoding='utf-8')
except AttributeError:
    pass

load_dotenv(Path("C:/Users/jmj2z/Projects/intelligence/ideas-hub/.env"))
addr = os.getenv("GMAIL_ADDRESS", "")
pwd = os.getenv("GMAIL_APP_PASSWORD", "")
if not addr or not pwd:
    print("Missing GMAIL_ADDRESS or GMAIL_APP_PASSWORD"); sys.exit(1)

HTML_FILE = Path(__file__).resolve().parent.parent / "results" / "N392_SELECTIVITY_SCREEN_RESULTS.html"

# Read or build HTML
if HTML_FILE.exists():
    html = HTML_FILE.read_text(encoding="utf-8")
else:
    # Build inline from markdown report
    md_file = HTML_FILE.with_suffix(".md")
    print(f"HTML not found, using inline version")
    html = INLINE_HTML  # fallback

msg = MIMEMultipart("alternative")
msg["Subject"] = "N392 Selectivity Mutation Screen \u2014 Results (2026-03-26)"
msg["From"] = addr
msg["To"] = addr

msg.attach(MIMEText("See HTML version.", "plain"))
msg.attach(MIMEText(html, "html"))

print(f"Sending to {addr}...")
with smtplib.SMTP_SSL("smtp.gmail.com", 465) as server:
    server.login(addr, pwd)
    server.sendmail(addr, addr, msg.as_string())
print("Sent!")
