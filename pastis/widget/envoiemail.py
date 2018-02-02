
import smtplib

sender = 'noreply@oamp.fr'
receivers = 'hebrard@iap.fr'

message = """From: PASTIS is fuuuuun
Subject: PASTIS simulations status

The PASTIS run is successfully terminated.
Thank you to use PASTIS - use with(out) moderation


This email is generated automatically by PASTIS, please do not reply."""

try :
   smtpObj = smtplib.SMTP('smtps.oamp.fr')
   smtpObj.sendmail(sender, receivers, message)
   print "Successfully sent email"
except SMTPException:
   print "Error: unable to send email"
