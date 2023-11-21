FROM rocker/r-ubuntu:jammy

RUN apt-get update && apt-get -y install python3-pip

COPY ./install_packages.r /home/

COPY ./requirements.txt /home/

RUN Rscript /home/install_packages.r

RUN pip3 install -r /home/requirements.txt

RUN wget https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.1-linux-kernel-3-x86_64.zip && unzip gcta-1.94.1-linux-kernel-3-x86_64.zip && cp /gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 /usr/local/bin/gcta && rm *.zip 
