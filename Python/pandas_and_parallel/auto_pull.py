import datetime as dt
import subprocess, multiprocessing
import os

# To find list of files from within /data folder on LMA box
# Yesterday
# d = dt.date.today()
# d = d-dt.timedelta(days =1)
# year   = d.year
# month  = d.month
# day    = d.day

year   = 2012
month  = 05
day    = 03

hour   = 00
minute = 00 
length = 24   # Length of time needed (in hours)

length = length * 6
original = dt.datetime(year, month, day, hour, minute)

station_ids = [('ip', 'Newdl', 'N', '216.48.255.140'),
               ('ip', 'Reese', 'O', '129.118.100.135'), 
               ('ip', 'Roose', 'R', '216.48.252.19'),
               ('ip', 'Peter', 'P', '216.48.250.178'), 
               ('ip', 'Loren', 'L', '216.48.249.100'),
               # ('ip', 'Idalo', 'G', '166.148.147.182'),
               # ('ip', 'Llano', 'W', '166.239.230.173'),
               # ('port', 'Biggi', 'B', '17095'), 
               ('port', 'Estac', 'E', '17070')
               ]


def new_list(name, symbol):
	list = ["nonsense"]
	list.append('/data/log/T%s%s' %(symbol, original.strftime('%y%m%d')))
	for i in range(length):
		new = original + i*(dt.timedelta(minutes=10))
		list.append('/data/%s/L%s_WestTexas_%s_%s_%s00.dat' % (new.strftime('%y%m%d'), symbol, name, new.strftime('%y%m%d'), new.strftime('%H%M')))
	list.append('more_nonsense')
	b = ' '.join(list)
	return b

def scp_ips(tuple_list0, tuple_list1, tuple_list2, tuple_list3):
	b = new_list(tuple_list1, tuple_list2)
	if tuple_list0 == 'ip':
		subprocess.call(['scp', 'root@%s:%s' %(tuple_list3, '\\\"'+b+'\\\"'), '/data/rtlma/exchange/raw/20%s' %(original.strftime('%y/%b/%d'))])
	if tuple_list0 == 'port':	
		subprocess.call(['sudo', 'scp', '-P', tuple_list3, 'root@localhost:%s' %('\\\"'+b+'\\\"'), '/data/rtlma/exchange/raw/20%s' %(original.strftime('%y/%b/%d'))])

jobs = []
if __name__ == '__main__':
	
	base_out_dir = ("/data/rtlma/exchange/raw/20%s" %(original.strftime('%y/%b/%d')))
	if os.path.exists(base_out_dir) == False:
		os.makedirs(base_out_dir)
		subprocess.call(['chmod', 'a+w', base_out_dir, '/data/rtlma/exchange/raw/20%s' %(date.strftime('%y/%b')), '/data/rtlma/exchange/raw/20%s' %(date.strftime('%y'))])

	for item in station_ids:
		proc = multiprocessing.Process(target = scp_ips, args = (item))
		jobs.append(proc)
		proc.start()