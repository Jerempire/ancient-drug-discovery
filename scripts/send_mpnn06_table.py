"""Send mpnn06 results table via Gmail."""
import sys, smtplib, os
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from dotenv import load_dotenv
from pathlib import Path

try: sys.stdout.reconfigure(encoding='utf-8')
except: pass

load_dotenv(Path('C:/Users/jmj2z/Projects/intelligence/substack-digest/.env'))
pw = os.environ.get('GMAIL_APP_PASSWORD', '')
recipient = 'jerempire99@gmail.com'

html = """
<div style="font-family: Georgia, serif; max-width: 800px; margin: 0 auto;">
<h1 style="border-bottom: 3px solid #2c5f2d; padding-bottom: 8px;">mpnn06 — Full Results Summary</h1>

<h2>Lead Construct: mpnn06</h2>
<table border="1" cellpadding="6" cellspacing="0" style="border-collapse:collapse; font-size:13px;">
<tr style="background:#2c5f2d; color:white;"><th>Property</th><th>Value</th></tr>
<tr><td>Size</td><td>113aa, ~12.4 kDa, single chain</td></tr>
<tr><td>ERAP2 ipTM</td><td><strong>0.815</strong> [0.878, 0.821, 0.745]</td></tr>
<tr><td>ERAP1 ipTM</td><td>0.326 (<strong>2.5x selectivity</strong>)</td></tr>
<tr><td>IRAP ipTM</td><td>0.454 (<strong>1.8x selectivity</strong>) [0.710, 0.346, 0.305]</td></tr>
<tr><td>Mutations (14)</td><td>R3I Y5H F6Q W19K L23M T37N I39T V41L F51D I56Y A89N L90A K91L N92R</td></tr>
<tr><td>Architecture</td><td>VAGSAF &mdash; (EAAAA)x3 helix &mdash; redesigned binder body</td></tr>
</table>

<h2>Full Campaign Results (All Constructs)</h2>
<table border="1" cellpadding="5" cellspacing="0" style="border-collapse:collapse; font-size:12px; width:100%;">
<tr style="background:#f0f0f0;">
<th>Round</th><th>Construct</th><th>ERAP2 avg</th><th>ERAP2 samples</th><th>ERAP1 avg</th><th>IRAP avg</th><th>IRAP samples</th><th>E2/E1</th><th>E2/IR</th>
</tr>
<tr><td>R1</td><td>Parent (92aa)</td><td>0.572</td><td>[0.650, 0.602, 0.463]</td><td>0.604</td><td>0.407</td><td>[0.465, 0.405, 0.350]</td><td style="color:red;">0.9x</td><td>1.4x</td></tr>
<tr><td>R1</td><td>B (V19W+I49W)</td><td>0.507</td><td>[0.661, 0.454, 0.406]</td><td>0.623</td><td>0.370</td><td>[0.386, 0.366, 0.357]</td><td style="color:red;">0.8x</td><td>1.4x</td></tr>
<tr><td>R1</td><td>C (GGGGS cargo)</td><td>0.507</td><td>[0.618, 0.469, 0.435]</td><td>0.518</td><td>0.424</td><td>[0.457, 0.433, 0.382]</td><td style="color:red;">1.0x</td><td>1.2x</td></tr>
<tr><td>R1</td><td>D (combined)</td><td>0.455</td><td>[0.522, 0.424, 0.418]</td><td>0.616</td><td>0.376</td><td>[0.402, 0.377, 0.349]</td><td style="color:red;">0.7x</td><td>1.2x</td></tr>
<tr style="background:#d4edda;"><td>R3</td><td><strong>C2 EAAAK</strong></td><td><strong>0.759</strong></td><td>[0.807, 0.810, 0.659]</td><td><strong>0.193</strong></td><td>0.609</td><td>[0.709, 0.597, 0.520]</td><td style="color:green;"><strong>3.9x</strong></td><td>1.2x</td></tr>
<tr><td>R3</td><td>VAGSAF ctrl</td><td>0.763</td><td>[0.846, 0.733, 0.710]</td><td>0.565</td><td>0.782</td><td>[0.820, 0.766, 0.760]</td><td>1.4x</td><td style="color:red;">1.0x</td></tr>
<tr><td>R3</td><td>VAGSAFLL</td><td>0.817</td><td>[0.828, 0.819, 0.805]</td><td>0.599</td><td>0.804</td><td>[0.870, 0.795, 0.747]</td><td>1.4x</td><td style="color:red;">1.0x</td></tr>
<tr><td>R3</td><td>CVAGSAFC disulfide</td><td>0.772</td><td>[0.837, 0.789, 0.688]</td><td>0.755</td><td>0.811</td><td>[0.842, 0.824, 0.767]</td><td style="color:red;">1.0x</td><td style="color:red;">1.0x</td></tr>
<tr><td>R3</td><td>PVAGSAFP Pro-cap</td><td>0.742</td><td>[0.767, 0.759, 0.700]</td><td>0.589</td><td>0.806</td><td>[0.822, 0.794, 0.803]</td><td>1.3x</td><td style="color:red;">0.9x</td></tr>
<tr><td>R3</td><td>GVAGSAFG cyclic</td><td>0.711</td><td>[0.793, 0.720, 0.620]</td><td>0.502</td><td>0.897</td><td>[0.939, 0.865, 0.887]</td><td>1.4x</td><td style="color:red;">0.8x</td></tr>
<tr style="background:#d4edda;"><td>R4</td><td><strong>c2_eaaaa3</strong></td><td><strong>0.778</strong></td><td>[0.827, 0.798, 0.711]</td><td>0.333</td><td><strong>0.573</strong></td><td>[0.668, 0.557, 0.494]</td><td>2.3x</td><td><strong>1.4x</strong></td></tr>
<tr><td>R4</td><td>c2_eaaak3 ctrl</td><td>0.729</td><td>[0.761, 0.755, 0.672]</td><td>0.294</td><td>0.588</td><td>[0.687, 0.571, 0.507]</td><td>2.5x</td><td>1.2x</td></tr>
<tr><td>R4</td><td>c2_eaaak2</td><td>0.703</td><td>[0.841, 0.666, 0.602]</td><td>0.240</td><td>0.571</td><td>[0.691, 0.702, 0.319]</td><td>2.9x</td><td>1.2x</td></tr>
<tr><td>R4</td><td>c2_eaaak1 shortest</td><td>0.743</td><td>[0.778, 0.762, 0.688]</td><td>0.299</td><td>0.740</td><td>[0.863, 0.668, 0.690]</td><td>2.5x</td><td style="color:red;">1.0x</td></tr>
<tr><td>R4</td><td>c2_eaaaa2</td><td>0.575</td><td>[0.600, 0.664, 0.462]</td><td>0.292</td><td>0.691</td><td>[0.686, 0.754, 0.634]</td><td>2.0x</td><td style="color:red;">0.8x</td></tr>
<tr><td>R4b</td><td>neg1 (I2E)</td><td>0.664</td><td>[0.815, 0.628, 0.550]</td><td>0.255</td><td>0.708</td><td>[0.902, 0.580, 0.642]</td><td>2.6x</td><td style="color:red;">0.9x</td></tr>
<tr><td>R4b</td><td>neg2 (D1R+I2E)</td><td>0.753</td><td>[0.797, 0.754, 0.708]</td><td>0.201</td><td>0.782</td><td>[0.824, 0.747, 0.774]</td><td>3.7x</td><td style="color:red;">1.0x</td></tr>
<tr><td>R4b</td><td>neg3 (all 3)</td><td>0.739</td><td>[0.757, 0.716, 0.743]</td><td>0.234</td><td>0.673</td><td>[0.682, 0.693, 0.643]</td><td>3.2x</td><td>1.1x</td></tr>
<tr style="background:#c8e6c9;"><td><strong>R5</strong></td><td><strong>mpnn06</strong></td><td><strong>0.815</strong></td><td><strong>[0.878, 0.821, 0.745]</strong></td><td><strong>0.326</strong></td><td><strong>0.454</strong></td><td><strong>[0.710, 0.346, 0.305]</strong></td><td><strong>2.5x</strong></td><td style="color:green;"><strong>1.8x</strong></td></tr>
<tr><td>R5</td><td>mpnn01</td><td>0.810</td><td>[0.913, 0.797, 0.721]</td><td>0.381</td><td>0.711</td><td>[0.782, 0.719, 0.632]</td><td>2.1x</td><td>1.1x</td></tr>
<tr><td>R5</td><td>mpnn02</td><td>0.773</td><td>[0.818, 0.779, 0.721]</td><td>0.289</td><td>0.785</td><td>[0.767, 0.804, 0.784]</td><td>2.7x</td><td style="color:red;">1.0x</td></tr>
<tr><td>R5</td><td>mpnn03</td><td>0.727</td><td>[0.830, 0.723, 0.628]</td><td>0.567</td><td>0.824</td><td>[0.879, 0.782, 0.811]</td><td style="color:red;">1.3x</td><td style="color:red;">0.9x</td></tr>
<tr><td>R5</td><td>mpnn04</td><td>0.753</td><td>[0.807, 0.747, 0.707]</td><td>0.242</td><td>0.705</td><td>[0.730, 0.782, 0.603]</td><td>3.1x</td><td>1.1x</td></tr>
<tr><td>R5</td><td>mpnn05</td><td>0.730</td><td>[0.759, 0.736, 0.696]</td><td>0.454</td><td>0.671</td><td>[0.782, 0.695, 0.536]</td><td>1.6x</td><td>1.1x</td></tr>
<tr><td>R5</td><td>mpnn07</td><td>0.676</td><td>[0.797, 0.672, 0.561]</td><td>0.267</td><td>0.708</td><td>[0.717, 0.738, 0.671]</td><td>2.5x</td><td style="color:red;">1.0x</td></tr>
<tr><td>R5</td><td>mpnn08</td><td>0.668</td><td>[0.753, 0.645, 0.608]</td><td>0.300</td><td>0.706</td><td>[0.721, 0.722, 0.676]</td><td>2.2x</td><td style="color:red;">0.9x</td></tr>
<tr><td>R5</td><td>mpnn09</td><td>0.730</td><td>[0.818, 0.715, 0.656]</td><td>0.322</td><td>0.745</td><td>[0.823, 0.754, 0.658]</td><td>2.3x</td><td style="color:red;">1.0x</td></tr>
<tr><td>R5</td><td>mpnn10</td><td>0.765</td><td>[0.870, 0.729, 0.698]</td><td>0.400</td><td>0.757</td><td>[0.794, 0.757, 0.721]</td><td>1.9x</td><td style="color:red;">1.0x</td></tr>
</table>

<h2>Next: Tweak mpnn06 to Drop IRAP Below 0.30</h2>
<p>IRAP samples [0.710, 0.346, 0.305] are bimodal &mdash; 2 of 3 already at target. Plan:</p>
<ol>
<li>Regenerate CIFs + analyze the 0.710 outlier IRAP pose</li>
<li>Test shorter linker (EAAAA)x2 on mpnn06</li>
<li>Test cargo variants: VAGSAE (F&rarr;E), VAGSAK (F&rarr;K) at P1</li>
<li>Targeted mutations from contact analysis</li>
</ol>

<p style="color:#888;font-size:12px;margin-top:20px;">162 predictions | 5 rounds | ~$2 compute | Jerempire/ancient-drug-discovery</p>
</div>
"""

msg = MIMEMultipart('alternative')
msg['From'] = recipient
msg['To'] = recipient
msg['Subject'] = 'mpnn06 Full Results Table + Next Steps'
msg.attach(MIMEText('See HTML version for full table.', 'plain', 'utf-8'))
msg.attach(MIMEText(html, 'html', 'utf-8'))

with smtplib.SMTP('smtp.gmail.com', 587) as server:
    server.starttls()
    server.login(recipient, pw)
    server.sendmail(recipient, recipient, msg.as_string())
    print('Email sent successfully')
