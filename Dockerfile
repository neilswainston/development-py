FROM ubuntu:latest
RUN apt-get update -y
RUN apt-get install -y python-pip python-dev build-essential
COPY . /
WORKDIR /
ENTRYPOINT ["python"]
CMD ["test.py", "https://ice.synbiochem.co.uk/login"]