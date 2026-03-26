"""Send contingency plan report via Gmail SMTP."""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")
except AttributeError:
    pass

import os, smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from pathlib import Path
from dotenv import load_dotenv

load_dotenv(Path("C:/Users/jmj2z/Projects/intelligence/ideas-hub/.env"))
GMAIL = os.getenv("GMAIL_ADDRESS", "")
PASS = os.getenv("GMAIL_APP_PASSWORD", "")

SUBJECT = "ERAP2: Contingency Plan — What To Do If 6-mers Drift in MD"

HTML = """\
<!DOCTYPE html>
<html><body style="font-family: Georgia, serif; max-width: 760px; margin: 0 auto; line-height: 1.7; color: #222;">

<h1 style="border-bottom: 2px solid #c53030;">Contingency Plan: If 6-mer Peptides Drift in MD</h1>
<p style="color: #666; font-size: 14px;">March 25, 2026 — Prepared while MD runs in progress</p>

<h2>Current Status</h2>
<ul>
<li><strong>Gate PASSED:</strong> FASGAV (scrambled control) ejected from ERAP2-K392 — 28.9 A COM drift, zero contacts. Pocket is NOT sticky.</li>
<li><strong>VAGSAF drifted</strong> — but likely due to 6.06 A grafting error (bad starting placement), not genuine non-binding.</li>
<li><strong>IAFSAF running now</strong> — 1.16 A grafting RMSD (good placement). Result in ~1 hour. This is the definitive test.</li>
</ul>

<h2>If IAFSAF Also Drifts — The 6 Options (Priority Order)</h2>

<h3 style="color: #2563eb;">1. Restrained Equilibration (TRY FIRST)</h3>
<p><strong>Cost: ~$0.50 | Time: ~3 hrs | Probability: HIGH</strong></p>
<p>Current protocol has no initial restraints. Water jostle can eject a small 6-mer before it settles. Fix: hold peptide in place for 1-2 ns while water/protein relax, then gradually release. Standard practice for small ligand MD — we skipped it for speed.</p>
<p><strong>If this fixes it:</strong> binding is real, we just needed proper equilibration.</p>

<h3 style="color: #2563eb;">2. Fall Back to 10-11 mers</h3>
<p><strong>Cost: ~$2 | Time: ~8 hrs | Probability: HIGH</strong></p>
<p>KILKLYSSKKY (11-mer) scored 0.892 ipTM, +0.427 selectivity, KKY motif confirmed load-bearing. More surface area = more grip. Trade-off: lose ERAP1 molecular ruler evasion, need sequence-based ERAP1 selectivity instead.</p>

<h3 style="color: #2563eb;">3. Cyclize the 6-mer</h3>
<p><strong>Cost: ~$500-800 synthesis | Time: 1-2 weeks | Probability: MEDIUM</strong></p>
<p>Linear 6-mer is floppy — high entropy penalty. Head-to-tail cyclic peptide is pre-organized. Better binding free energy, more drug-like, protease-resistant. But Boltz-2 may not handle cyclic peptides well.</p>

<h3 style="color: #6b7280;">4. D-amino Acid Backbone</h3>
<p><strong>Cost: ~$200 | Time: ~1 week | Probability: MEDIUM</strong></p>
<p>Flip backbone geometry. Protease-resistant. Different conformations may fit better. Can model computationally first with PepINVENT.</p>

<h3 style="color: #6b7280;">5. Non-Natural Anchor Residue</h3>
<p><strong>Cost: ~$1,000+ | Time: 2-4 weeks | Probability: MEDIUM-LOW</strong></p>
<p>Boronic acid at P1 for reversible covalent zinc binding. N-methyl amino acids for reduced flexibility. Phase 2 optimization only.</p>

<h3 style="color: #6b7280;">6. Re-run All 18 with Restrained Equilibration</h3>
<p><strong>Cost: ~$6 | Time: ~36 hrs | Probability: MEDIUM</strong></p>
<p>Before giving up on 6-mers entirely, re-run full panel with proper protocol. Some may be stable even if IAFSAF isn't.</p>

<h2>Decision Tree</h2>
<pre style="background: #f8f8f8; padding: 15px; border-radius: 5px; font-size: 13px;">
IAFSAF drifts in MD
    |
    +--&gt; Try restrained equilibration (Option 1, ~$0.50)
    |       |
    |       +--&gt; Stays --&gt; 6-mers work! Protocol was the issue.
    |       |
    |       +--&gt; Still drifts --&gt; 6-mers don't bind stably
    |               |
    |               +--&gt; Run KILKLYSSKKY 11-mer (Option 2, ~$2)
    |               |       |
    |               |       +--&gt; Stays --&gt; Pivot to 11-mers
    |               |       +--&gt; Drifts --&gt; Deeper pose problem
    |               |
    |               +--&gt; Cyclize IAFSAF (Option 3, ~$800)
    |                       if 6-mer length is critical for ERAP1 evasion
    |
    +--&gt; Total cost to reach answer: $0.50 - $3.00
</pre>

<h2>What Happens If IAFSAF STAYS</h2>
<p>Then we don't need this plan. IAFSAF becomes the confirmed lead — a 6-mer that:</p>
<ul>
<li>Binds K392-ERAP2 stably (MD confirmed)</li>
<li>Evades ERAP1 (molecular ruler, Boltz-2 confirmed)</li>
<li>Evades IRAP (surface area evasion, Boltz-2 confirmed)</li>
<li>Ejects scrambled control (FASGAV gate passed)</li>
</ul>
<p>That's the computational case for synthesis ($100-150) and patent filing ($75).</p>

<h2>Key Lesson from VAGSAF</h2>
<p>VAGSAF drifted not because it doesn't bind, but because the grafting placed it 6 A from the correct position. <strong>Starting pose quality determines MD outcome.</strong> IAFSAF's 1.16 A graft alignment is much better. If it stays, we know the binding is real AND that careful structure preparation is essential. VAGSAF just needs a better starting pose (re-dock on full protein or use restrained equilibration).</p>

<p style="color: #666; font-size: 12px; margin-top: 40px; border-top: 1px solid #ccc; padding-top: 10px;">Ancient Drug Discovery Pipeline — March 25, 2026<br>Contingency plan saved at docs/CONTINGENCY_IF_6MER_DRIFTS.md</p>
</body></html>"""

def main():
    if not GMAIL or not PASS:
        print("ERROR: Gmail credentials not found"); sys.exit(1)
    msg = MIMEMultipart("alternative")
    msg["Subject"] = SUBJECT
    msg["From"] = GMAIL
    msg["To"] = GMAIL
    msg.attach(MIMEText(HTML, "html", "utf-8"))
    with smtplib.SMTP_SSL("smtp.gmail.com", 465) as s:
        s.login(GMAIL, PASS)
        s.sendmail(GMAIL, [GMAIL], msg.as_string())
    print("Email sent to %s" % GMAIL)

if __name__ == "__main__":
    main()
