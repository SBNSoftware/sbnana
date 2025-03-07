# Import smtplib for the actual sending function
import smtplib
# Import the email modules we'll need
from email.mime.text import MIMEText

import random
import sys
import shutil


files = sys.argv[1:]
random.shuffle(files)
OUTDIR = "/icarus/data/users/gputnam/DMCP2023G/mc-F-spectra/blinded-datas/"
alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
me = "gray@putnam"
you = "joanli411@gmail.com"

for i, f in enumerate(files):
    fname = f.split("/")[-1]
    src = f
    dst = OUTDIR + alphabet[i] + ".df"
    shutil.copyfile(src, dst)

    if fname == "thedata.df": # the data
        msg = MIMEText("The real data is in dataset %s." % alphabet[i])
        msg['Subject'] = 'This is the blinding email. Dont tell anyone!'
        msg['From'] = me
        msg['To'] = you

        # Send the message via our own SMTP server, but don't include the
        # envelope header.
        s = smtplib.SMTP('localhost')
        s.sendmail(me, [you], msg.as_string())
        s.quit()
