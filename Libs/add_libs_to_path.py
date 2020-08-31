# site_packages_path is the packages folder, which in my case is:
site_packages_path = r'C:\Users\Dhwani\AppData\Local\Continuum\anaconda3\Lib\site-packages'
site_packages_path = r'C:\Users\Hermínio\Appdata\Local\Programs\Python\Python38-32\Lib\site-packages'


# path that you wanna add, which again in my case is 
path_to_add = "C:\cygwin64\home\Libs\GeoFlow1D\Libs"

f = open(site_packages_path + "\custom_path.pth", "a")
f.write(path_to_add)