# Import smtplib for the actual sending function
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText


def email_status(user_email, run_name, run_id, message):
    # our_address = "recblast@reciprocalblast.com"
    # our_address = "moranneuhof@mail.tau.ac.il"

    # TODO: make this more secure
    googleapps_user = 'recblast@gmail.com'
    googleapps_pass = 'Recipr0ca1y'

    BCC = ['efratefrat120@gmail.com', 'neuhofmo@gmail.com']
    msg = MIMEMultipart()
    msg['Subject'] = "Your job {0} (run id: {1}) status update".format(run_name, run_id)
    msg['From'] = googleapps_user
    msg['To'] = user_email

    with open("templates/email_template.html", 'r') as f:
        html = f.read()
    # print(html)
    msg.attach(MIMEText(html.format(message), 'html'))

    server = smtplib.SMTP_SSL('smtp.gmail.com', 465)
    server.login(googleapps_user, googleapps_pass)
    server.sendmail(googleapps_user, [user_email] + BCC, msg.as_string())
    server.quit()
    return True
