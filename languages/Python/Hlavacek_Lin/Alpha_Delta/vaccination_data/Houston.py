from bs4 import BeautifulSoup 
from selenium import webdriver
import re
import numpy as np
url = [0]*9
y = np.zeros((len(url),650))
S0 = 7066141
url[0] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/austin-county/48015/'
url[1] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/brazoria-county/48039/'
url[2] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/chambers-county/48071/'
url[3] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/fortbend-county/48157/'
url[4] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/galveston-county/48167/'
url[5] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/harris-county/48201/'
url[6] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/liberty-county/48291/'
url[7] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/montgomery-county/48339/'
url[8] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/waller-county/48473/'
for ii in range(len(url)):
    driver = webdriver.Chrome()
    driver.get(url[ii])
    delay = 10
    soup = BeautifulSoup(driver.page_source, 'html.parser')
    driver.quit()
    with open('h'+str(ii)+'.txt', 'w') as f:
        f.write(str(soup))
    text_file = open('h'+str(ii)+'.txt','r')
    file_lines = text_file.readlines()
    temp = file_lines[1538][10:-2]
    exec(temp)
    x = []
    for jj in range(len(bcd_data)):
        x.append(bcd_data[jj]['cv'])
    temp1 = np.asarray(x)
    y[ii,350:] = temp1[:300]
    
temp2 = np.load('vac.npy')
temp2[4,:] = np.diff(np.sum(y,0))/S0
np.save('vac.npy',temp2)