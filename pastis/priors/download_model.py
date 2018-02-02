#! /usr/bin/python2.7
def besancon_model_form(name_target, model, minimum_distance, maximum_distance, longit, latit, solid_angle, maximum_magnitude, minimum_magnitude, email='', outfile = '/data/pasi/PASTIS/datafiles/GalModel.dat'):
    if model == 'BES':
        """
        Function that fills in the forms for the Besancon model with our data, more parameters can be added by looking at the source code of the web page.
        For the parameters that are not defined the model provides default values wich can be seen on the website of the model
        """
        import os
        os.chdir(os.path.split(outfile)[0])
        import mechanize
        br = mechanize.Browser()      # open the module browser
        br.set_handle_refresh(False)
        br.addheaders = [('User-agent', 'Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.9.0.1) Gecko/2008071615 Fedora/3.0.1-1.fc9 Firefox/3.0.1')]
        br.set_handle_robots(False)  # no robots
        br.open("http://model.obs-besancon.fr/modele_form.php?kleb=3&klec=1&cinem=0")  # connection to the model site
        br.select_form(nr=0)  # selection of the first form in the model
        br.form["email"] = email
        br.form["rinf"] = str(minimum_distance)
        br.form["rsup"] = str(maximum_distance)
        br.form["longit"] = str(longit)
        br.form["latit"] = str(latit)
        br.form["soli"] = str(solid_angle)
        br.form["band0[0]"] = "-99.00"
        br.form["bandf[0]"] = "99.00"
        br.form["band0[5]"] = str(minimum_magnitude)
        br.form["bandf[5]"] = str(maximum_magnitude)
        submit_response = br.submit(VALUE='submit')
        html_source = br.response().read()  # we save and return the web page source after the form has been submited wich gives us our file number
        print(html_source)
        from bs4 import BeautifulSoup as bs  # allows the extraction of data from web pages
        import re  # module that allows to work on text and extract the useful data

        soup = bs(html_source)  # saves the source as a string variable
        span = soup.find_all('span')  # we take all the text in the /span parts
        span_str = str(span[1])
        search_num = re.search('enhanced">(.+?)</span>', span_str)  # the file number is between the enhanced"> and the </span>, this will work as long as the besancon model html source is not changed by the besancon team
        print search_num
        if search_num:
            file_number = search_num.group(1)  # Not sure if this if is necessary !!!!

        import time
        import ftplib

        """
        Function that downloads the file containing the data
        """
        test2 = True
        server = 'sasftp.obs-besancon.fr'
        username = 'anonymous'
        directory = '/modele/modele2003'
        ftp = ftplib.FTP(server)
        ftp.login(username, email)
        ftp.cwd(directory)
        print("Connected")
        fhandle = open(file_number, 'wb')
        print("File created")
        while test2:
            try:
                print("Trying to download file in 30 seconds")
                time.sleep(30)
                print("Trying to download now")
                ftp.retrbinary('RETR ' + file_number, fhandle.write)
                fhandle.close()
                ftp.quit()
                test2 = False
                print("succees")
            except:
                test2 = True
                print("erreur")
        import os
        new_file_name = os.path.split(outfile)[1]
        for filename in os.listdir("."):
            if filename.startswith(file_number):
                os.rename(filename, new_file_name)
        return 
