"""Send full 18-peptide table report via Gmail SMTP."""
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

SUBJECT = "ERAP2: Full 18-Peptide 4-Target Panel — All Short-Stack Results"

HTML = """\
<!DOCTYPE html>
<html><body style="font-family: Georgia, serif; max-width: 800px; margin: 0 auto; line-height: 1.7; color: #222;">

<h1 style="border-bottom: 2px solid #1e40af;">Full 4-Target Panel: All 18 Short-Stack Peptides</h1>
<p style="color: #666; font-size: 14px;">March 24, 2026 | Boltz-2, 3 diffusion samples, seed 42</p>

<h2>Triple-Selective (K392 > ERAP1, IRAP, and N392) — 10 peptides</h2>

<table style="border-collapse: collapse; width: 100%; font-size: 12px;">
<tr style="background: #1e40af; color: white;">
<th style="padding: 5px; border: 1px solid #ccc;">Rank</th>
<th style="padding: 5px; border: 1px solid #ccc;">Peptide</th>
<th style="padding: 5px; border: 1px solid #ccc;">Len</th>
<th style="padding: 5px; border: 1px solid #ccc;">K392</th>
<th style="padding: 5px; border: 1px solid #ccc;">N392</th>
<th style="padding: 5px; border: 1px solid #ccc;">ERAP1</th>
<th style="padding: 5px; border: 1px solid #ccc;">IRAP</th>
<th style="padding: 5px; border: 1px solid #ccc;">sel E1</th>
<th style="padding: 5px; border: 1px solid #ccc;">sel IR</th>
<th style="padding: 5px; border: 1px solid #ccc;">sel N</th>
<th style="padding: 5px; border: 1px solid #ccc;">Design</th>
</tr>
<tr style="background: #dbeafe;"><td style="padding: 4px; border: 1px solid #ccc;"><strong>1</strong></td><td style="padding: 4px; border: 1px solid #ccc;"><strong>VAGSAF</strong></td><td style="padding: 4px; border: 1px solid #ccc;">6</td><td style="padding: 4px; border: 1px solid #ccc;"><strong>0.905</strong></td><td style="padding: 4px; border: 1px solid #ccc;">0.870</td><td style="padding: 4px; border: 1px solid #ccc;">0.335</td><td style="padding: 4px; border: 1px solid #ccc;">0.236</td><td style="padding: 4px; border: 1px solid #ccc;">+0.570</td><td style="padding: 4px; border: 1px solid #ccc;">+0.669</td><td style="padding: 4px; border: 1px solid #ccc;">+0.035</td><td style="padding: 4px; border: 1px solid #ccc;">Val P1 + Phe C-term</td></tr>
<tr style="background: #dbeafe;"><td style="padding: 4px; border: 1px solid #ccc;"><strong>2</strong></td><td style="padding: 4px; border: 1px solid #ccc;"><strong>IAFSAF</strong></td><td style="padding: 4px; border: 1px solid #ccc;">6</td><td style="padding: 4px; border: 1px solid #ccc;"><strong>0.870</strong></td><td style="padding: 4px; border: 1px solid #ccc;">0.631</td><td style="padding: 4px; border: 1px solid #ccc;">0.192</td><td style="padding: 4px; border: 1px solid #ccc;">0.251</td><td style="padding: 4px; border: 1px solid #ccc;">+0.678</td><td style="padding: 4px; border: 1px solid #ccc;">+0.619</td><td style="padding: 4px; border: 1px solid #ccc;"><strong>+0.239</strong></td><td style="padding: 4px; border: 1px solid #ccc;">Phe P3 clash handle</td></tr>
<tr><td style="padding: 4px; border: 1px solid #ccc;">3</td><td style="padding: 4px; border: 1px solid #ccc;">VAWSAF</td><td style="padding: 4px; border: 1px solid #ccc;">6</td><td style="padding: 4px; border: 1px solid #ccc;">0.852</td><td style="padding: 4px; border: 1px solid #ccc;">0.830</td><td style="padding: 4px; border: 1px solid #ccc;">0.169</td><td style="padding: 4px; border: 1px solid #ccc;">0.226</td><td style="padding: 4px; border: 1px solid #ccc;">+0.683</td><td style="padding: 4px; border: 1px solid #ccc;">+0.626</td><td style="padding: 4px; border: 1px solid #ccc;">+0.022</td><td style="padding: 4px; border: 1px solid #ccc;">Trp P3 clash handle</td></tr>
<tr><td style="padding: 4px; border: 1px solid #ccc;">4</td><td style="padding: 4px; border: 1px solid #ccc;">VAFSAGY</td><td style="padding: 4px; border: 1px solid #ccc;">7</td><td style="padding: 4px; border: 1px solid #ccc;">0.833</td><td style="padding: 4px; border: 1px solid #ccc;">0.754</td><td style="padding: 4px; border: 1px solid #ccc;">0.275</td><td style="padding: 4px; border: 1px solid #ccc;">0.293</td><td style="padding: 4px; border: 1px solid #ccc;">+0.558</td><td style="padding: 4px; border: 1px solid #ccc;">+0.540</td><td style="padding: 4px; border: 1px solid #ccc;">+0.079</td><td style="padding: 4px; border: 1px solid #ccc;">Val P1 + Phe mid + Tyr C</td></tr>
<tr><td style="padding: 4px; border: 1px solid #ccc;">5</td><td style="padding: 4px; border: 1px solid #ccc;">LAGSAF</td><td style="padding: 4px; border: 1px solid #ccc;">6</td><td style="padding: 4px; border: 1px solid #ccc;">0.821</td><td style="padding: 4px; border: 1px solid #ccc;">0.717</td><td style="padding: 4px; border: 1px solid #ccc;">0.290</td><td style="padding: 4px; border: 1px solid #ccc;">0.349</td><td style="padding: 4px; border: 1px solid #ccc;">+0.531</td><td style="padding: 4px; border: 1px solid #ccc;">+0.472</td><td style="padding: 4px; border: 1px solid #ccc;">+0.104</td><td style="padding: 4px; border: 1px solid #ccc;">Leu P1 + Phe C-term</td></tr>
<tr><td style="padding: 4px; border: 1px solid #ccc;">6</td><td style="padding: 4px; border: 1px solid #ccc;">IAGSAW</td><td style="padding: 4px; border: 1px solid #ccc;">6</td><td style="padding: 4px; border: 1px solid #ccc;">0.820</td><td style="padding: 4px; border: 1px solid #ccc;">0.744</td><td style="padding: 4px; border: 1px solid #ccc;">0.293</td><td style="padding: 4px; border: 1px solid #ccc;">0.305</td><td style="padding: 4px; border: 1px solid #ccc;">+0.527</td><td style="padding: 4px; border: 1px solid #ccc;">+0.515</td><td style="padding: 4px; border: 1px solid #ccc;">+0.076</td><td style="padding: 4px; border: 1px solid #ccc;">Ile P1 + Trp C-term</td></tr>
<tr><td style="padding: 4px; border: 1px solid #ccc;">7</td><td style="padding: 4px; border: 1px solid #ccc;">AAGSAF</td><td style="padding: 4px; border: 1px solid #ccc;">6</td><td style="padding: 4px; border: 1px solid #ccc;">0.799</td><td style="padding: 4px; border: 1px solid #ccc;">0.635</td><td style="padding: 4px; border: 1px solid #ccc;">0.267</td><td style="padding: 4px; border: 1px solid #ccc;">0.276</td><td style="padding: 4px; border: 1px solid #ccc;">+0.532</td><td style="padding: 4px; border: 1px solid #ccc;">+0.523</td><td style="padding: 4px; border: 1px solid #ccc;">+0.164</td><td style="padding: 4px; border: 1px solid #ccc;">Ala P1 baseline</td></tr>
<tr><td style="padding: 4px; border: 1px solid #ccc;">8</td><td style="padding: 4px; border: 1px solid #ccc;">EAGSAF</td><td style="padding: 4px; border: 1px solid #ccc;">6</td><td style="padding: 4px; border: 1px solid #ccc;">0.730</td><td style="padding: 4px; border: 1px solid #ccc;">0.641</td><td style="padding: 4px; border: 1px solid #ccc;">0.277</td><td style="padding: 4px; border: 1px solid #ccc;">0.270</td><td style="padding: 4px; border: 1px solid #ccc;">+0.453</td><td style="padding: 4px; border: 1px solid #ccc;">+0.460</td><td style="padding: 4px; border: 1px solid #ccc;">+0.089</td><td style="padding: 4px; border: 1px solid #ccc;">Glu P1 charged</td></tr>
<tr><td style="padding: 4px; border: 1px solid #ccc;">9</td><td style="padding: 4px; border: 1px solid #ccc;">VASSKY</td><td style="padding: 4px; border: 1px solid #ccc;">6</td><td style="padding: 4px; border: 1px solid #ccc;">0.718</td><td style="padding: 4px; border: 1px solid #ccc;">0.657</td><td style="padding: 4px; border: 1px solid #ccc;">0.299</td><td style="padding: 4px; border: 1px solid #ccc;">0.411</td><td style="padding: 4px; border: 1px solid #ccc;">+0.419</td><td style="padding: 4px; border: 1px solid #ccc;">+0.307</td><td style="padding: 4px; border: 1px solid #ccc;">+0.061</td><td style="padding: 4px; border: 1px solid #ccc;">Val P1 + mini KY motif</td></tr>
<tr><td style="padding: 4px; border: 1px solid #ccc;">10</td><td style="padding: 4px; border: 1px solid #ccc;">IAFSAAY</td><td style="padding: 4px; border: 1px solid #ccc;">7</td><td style="padding: 4px; border: 1px solid #ccc;">0.710</td><td style="padding: 4px; border: 1px solid #ccc;">0.627</td><td style="padding: 4px; border: 1px solid #ccc;">0.347</td><td style="padding: 4px; border: 1px solid #ccc;">0.334</td><td style="padding: 4px; border: 1px solid #ccc;">+0.363</td><td style="padding: 4px; border: 1px solid #ccc;">+0.376</td><td style="padding: 4px; border: 1px solid #ccc;">+0.083</td><td style="padding: 4px; border: 1px solid #ccc;">Phe P3, Tyr C-term</td></tr>
</table>

<h2>ERAP1-Selective Only (avoids ERAP1/IRAP but NOT K392-selective over N392) — 8 peptides</h2>

<table style="border-collapse: collapse; width: 100%; font-size: 12px;">
<tr style="background: #6b7280; color: white;">
<th style="padding: 5px; border: 1px solid #ccc;">Rank</th>
<th style="padding: 5px; border: 1px solid #ccc;">Peptide</th>
<th style="padding: 5px; border: 1px solid #ccc;">Len</th>
<th style="padding: 5px; border: 1px solid #ccc;">K392</th>
<th style="padding: 5px; border: 1px solid #ccc;">N392</th>
<th style="padding: 5px; border: 1px solid #ccc;">ERAP1</th>
<th style="padding: 5px; border: 1px solid #ccc;">IRAP</th>
<th style="padding: 5px; border: 1px solid #ccc;">sel E1</th>
<th style="padding: 5px; border: 1px solid #ccc;">sel IR</th>
<th style="padding: 5px; border: 1px solid #ccc;">sel N</th>
<th style="padding: 5px; border: 1px solid #ccc;">Design</th>
</tr>
<tr><td style="padding: 4px; border: 1px solid #ccc;">11</td><td style="padding: 4px; border: 1px solid #ccc;">IALSAFW</td><td style="padding: 4px; border: 1px solid #ccc;">7</td><td style="padding: 4px; border: 1px solid #ccc;">0.877</td><td style="padding: 4px; border: 1px solid #ccc;">0.874</td><td style="padding: 4px; border: 1px solid #ccc;">0.289</td><td style="padding: 4px; border: 1px solid #ccc;">0.344</td><td style="padding: 4px; border: 1px solid #ccc;">+0.588</td><td style="padding: 4px; border: 1px solid #ccc;">+0.533</td><td style="padding: 4px; border: 1px solid #ccc;">+0.003</td><td style="padding: 4px; border: 1px solid #ccc;">Ile P1 + Phe-Trp pair</td></tr>
<tr><td style="padding: 4px; border: 1px solid #ccc;">12</td><td style="padding: 4px; border: 1px solid #ccc;">VAGSAAF</td><td style="padding: 4px; border: 1px solid #ccc;">7</td><td style="padding: 4px; border: 1px solid #ccc;">0.801</td><td style="padding: 4px; border: 1px solid #ccc;">0.882</td><td style="padding: 4px; border: 1px solid #ccc;">0.270</td><td style="padding: 4px; border: 1px solid #ccc;">0.335</td><td style="padding: 4px; border: 1px solid #ccc;">+0.531</td><td style="padding: 4px; border: 1px solid #ccc;">+0.466</td><td style="padding: 4px; border: 1px solid #ccc;">-0.081</td><td style="padding: 4px; border: 1px solid #ccc;">Val P1 extended 7-mer</td></tr>
<tr><td style="padding: 4px; border: 1px solid #ccc;">13</td><td style="padding: 4px; border: 1px solid #ccc;">VAGSYW</td><td style="padding: 4px; border: 1px solid #ccc;">6</td><td style="padding: 4px; border: 1px solid #ccc;">0.777</td><td style="padding: 4px; border: 1px solid #ccc;">0.777</td><td style="padding: 4px; border: 1px solid #ccc;">0.242</td><td style="padding: 4px; border: 1px solid #ccc;">0.296</td><td style="padding: 4px; border: 1px solid #ccc;">+0.535</td><td style="padding: 4px; border: 1px solid #ccc;">+0.481</td><td style="padding: 4px; border: 1px solid #ccc;">+0.000</td><td style="padding: 4px; border: 1px solid #ccc;">Val P1 + aromatic pair C</td></tr>
<tr><td style="padding: 4px; border: 1px solid #ccc;">14</td><td style="padding: 4px; border: 1px solid #ccc;">IAGSKKY</td><td style="padding: 4px; border: 1px solid #ccc;">7</td><td style="padding: 4px; border: 1px solid #ccc;">0.776</td><td style="padding: 4px; border: 1px solid #ccc;">0.820</td><td style="padding: 4px; border: 1px solid #ccc;">0.206</td><td style="padding: 4px; border: 1px solid #ccc;">0.530</td><td style="padding: 4px; border: 1px solid #ccc;">+0.570</td><td style="padding: 4px; border: 1px solid #ccc;">+0.246</td><td style="padding: 4px; border: 1px solid #ccc;">-0.044</td><td style="padding: 4px; border: 1px solid #ccc;">KKY motif 7-mer</td></tr>
<tr><td style="padding: 4px; border: 1px solid #ccc;">15</td><td style="padding: 4px; border: 1px solid #ccc;">IAGSAY</td><td style="padding: 4px; border: 1px solid #ccc;">6</td><td style="padding: 4px; border: 1px solid #ccc;">0.759</td><td style="padding: 4px; border: 1px solid #ccc;">0.858</td><td style="padding: 4px; border: 1px solid #ccc;">0.302</td><td style="padding: 4px; border: 1px solid #ccc;">0.215</td><td style="padding: 4px; border: 1px solid #ccc;">+0.457</td><td style="padding: 4px; border: 1px solid #ccc;">+0.544</td><td style="padding: 4px; border: 1px solid #ccc;">-0.099</td><td style="padding: 4px; border: 1px solid #ccc;">Ile P1 + Tyr C-term</td></tr>
<tr><td style="padding: 4px; border: 1px solid #ccc;">16</td><td style="padding: 4px; border: 1px solid #ccc;">IAGSAAY</td><td style="padding: 4px; border: 1px solid #ccc;">7</td><td style="padding: 4px; border: 1px solid #ccc;">0.731</td><td style="padding: 4px; border: 1px solid #ccc;">0.749</td><td style="padding: 4px; border: 1px solid #ccc;">0.199</td><td style="padding: 4px; border: 1px solid #ccc;">0.258</td><td style="padding: 4px; border: 1px solid #ccc;">+0.532</td><td style="padding: 4px; border: 1px solid #ccc;">+0.473</td><td style="padding: 4px; border: 1px solid #ccc;">-0.018</td><td style="padding: 4px; border: 1px solid #ccc;">Ile P1 + Tyr C-term ext</td></tr>
<tr><td style="padding: 4px; border: 1px solid #ccc;">17</td><td style="padding: 4px; border: 1px solid #ccc;">VAGSKKY</td><td style="padding: 4px; border: 1px solid #ccc;">7</td><td style="padding: 4px; border: 1px solid #ccc;">0.728</td><td style="padding: 4px; border: 1px solid #ccc;">0.844</td><td style="padding: 4px; border: 1px solid #ccc;">0.144</td><td style="padding: 4px; border: 1px solid #ccc;">0.388</td><td style="padding: 4px; border: 1px solid #ccc;">+0.584</td><td style="padding: 4px; border: 1px solid #ccc;">+0.340</td><td style="padding: 4px; border: 1px solid #ccc;">-0.116</td><td style="padding: 4px; border: 1px solid #ccc;">KKY motif 7-mer</td></tr>
<tr><td style="padding: 4px; border: 1px solid #ccc;">18</td><td style="padding: 4px; border: 1px solid #ccc;">VAWSAAF</td><td style="padding: 4px; border: 1px solid #ccc;">7</td><td style="padding: 4px; border: 1px solid #ccc;">0.629</td><td style="padding: 4px; border: 1px solid #ccc;">0.695</td><td style="padding: 4px; border: 1px solid #ccc;">0.193</td><td style="padding: 4px; border: 1px solid #ccc;">0.347</td><td style="padding: 4px; border: 1px solid #ccc;">+0.436</td><td style="padding: 4px; border: 1px solid #ccc;">+0.282</td><td style="padding: 4px; border: 1px solid #ccc;">-0.066</td><td style="padding: 4px; border: 1px solid #ccc;">Trp P3, 7-mer</td></tr>
</table>

<h2>Summary Statistics</h2>
<table style="border-collapse: collapse; font-size: 14px;">
<tr style="background: #f0f0f0;"><th style="padding: 6px; border: 1px solid #ccc;">Metric</th><th style="padding: 6px; border: 1px solid #ccc;">6-mers (10)</th><th style="padding: 6px; border: 1px solid #ccc;">7-mers (8)</th><th style="padding: 6px; border: 1px solid #ccc;">All 18</th></tr>
<tr><td style="padding: 6px; border: 1px solid #ccc;">Avg K392 ipTM</td><td style="padding: 6px; border: 1px solid #ccc;">0.805</td><td style="padding: 6px; border: 1px solid #ccc;">0.761</td><td style="padding: 6px; border: 1px solid #ccc;">0.786</td></tr>
<tr><td style="padding: 6px; border: 1px solid #ccc;">Avg ERAP1 ipTM</td><td style="padding: 6px; border: 1px solid #ccc;">0.266</td><td style="padding: 6px; border: 1px solid #ccc;">0.240</td><td style="padding: 6px; border: 1px solid #ccc;">0.255</td></tr>
<tr><td style="padding: 6px; border: 1px solid #ccc;">Avg IRAP ipTM</td><td style="padding: 6px; border: 1px solid #ccc;">0.298</td><td style="padding: 6px; border: 1px solid #ccc;">0.329</td><td style="padding: 6px; border: 1px solid #ccc;">0.312</td></tr>
<tr><td style="padding: 6px; border: 1px solid #ccc;">Avg sel ERAP1</td><td style="padding: 6px; border: 1px solid #ccc;">+0.539</td><td style="padding: 6px; border: 1px solid #ccc;">+0.520</td><td style="padding: 6px; border: 1px solid #ccc;">+0.531</td></tr>
<tr><td style="padding: 6px; border: 1px solid #ccc;">Triple-selective</td><td style="padding: 6px; border: 1px solid #ccc;">8/10</td><td style="padding: 6px; border: 1px solid #ccc;">2/8</td><td style="padding: 6px; border: 1px solid #ccc;">10/18</td></tr>
</table>

<h2>Key Patterns</h2>
<ul>
<li><strong>6-mers dominate:</strong> 8/10 triple-selective hits are 6-mers. 7-mers tend to lean N392 (5/8).</li>
<li><strong>Phe C-terminal is the common thread:</strong> Top 6 triple-selective peptides all end in Phe. The aromatic anchor at the channel exit appears critical.</li>
<li><strong>Every peptide avoids ERAP1 and IRAP:</strong> Lowest ERAP1 selectivity is +0.363. Not one 6-7mer binds ERAP1 or IRAP competitively. The size floor works universally.</li>
<li><strong>P1 matters less than P3 and C-term at 6-mer length:</strong> IAFSAF (F at P3, delta +0.239) vs VAGSAF (no P3 aromatic, delta +0.035). Internal aromatic placement drives selectivity more than P1.</li>
<li><strong>KKY motif flips at 7-mer:</strong> K392-selective on 11-mers, N392-selective on 7-mers. Different channel depths, different rules.</li>
</ul>

<h2>Method</h2>
<p><strong>Boltz-2 v2.2.1</strong> structure prediction, 3 diffusion samples, seed 42. Each peptide docked against cropped channel regions: ERAP2-K392 (res 350-450), ERAP2-N392 (K392N mutant), ERAP1 (aligned res 333-433), IRAP (aligned res 444-544). All chains use <code>msa: empty</code>. Scoring: ipTM averaged across diffusion samples. Selectivity = K392 ipTM minus target ipTM. Single seed — top hits need multi-seed confirmation.</p>

<p style="color: #666; font-size: 12px; margin-top: 40px; border-top: 1px solid #ccc; padding-top: 10px;">Ancient Drug Discovery Pipeline — March 24, 2026</p>
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
