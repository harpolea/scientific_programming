import socket
import silly
import time

def hello(s):
    print("Hello {}!".format(s))
    print(silly.sentence())

if __name__ == "__main__":
    name = socket.gethostname() # get system hostname
    while(True):
        hello(name)
        time.sleep(2)
