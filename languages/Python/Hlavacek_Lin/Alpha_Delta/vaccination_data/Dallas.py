from bs4 import BeautifulSoup 
from selenium import webdriver
import re
import numpy as np
url = [0]*11
y = np.zeros((len(url),650))
S0 = 7573136
url[0] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/collin-county/48085/'
url[1] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/dallas-county/48113/'
url[2] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/denton-county/48121/'
url[3] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/ellis-county/48139/'
url[4] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/hunt-county/48231/'
url[5] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/kaufman-county/48257/'
url[6] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/rockwall-county/48397/'
url[7] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/johnson-county/48251/'
url[8] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/parker-county/48367/'
url[9] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/tarrant-county/48439/'
url[10] = 'https://data.democratandchronicle.com/covid-19-vaccine-tracker/texas/wise-county/48497/'
for ii in range(len(url)):
    driver = webdriver.Chrome()
    driver.get(url[ii])
    delay = 10
    soup = BeautifulSoup(driver.page_source, 'html.parser')
    driver.quit()
    with open('d'+str(ii)+'.txt', 'w') as f:
        f.write(str(soup))
    text_file = open('d'+str(ii)+'.txt','r')
    file_lines = text_file.readlines()
    temp = file_lines[1538][10:-2]
    exec(temp)
    x = []
    for jj in range(len(bcd_data)):
        x.append(bcd_data[jj]['cv'])
    temp1 = np.asarray(x)
    y[ii,350:] = temp1[:300]
    
temp2 = np.load('vac.npy')
temp2[3,:] = np.diff(np.sum(y,0))/S0
np.save('vac.npy',temp2)