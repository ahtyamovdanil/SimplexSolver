# Get the base Ubuntu image from Docker Hub
FROM ubuntu:latest

# Update apps on the base image
RUN apt-get -y update && apt-get install -y

# Install mpich library
RUN apt-get -y install mpich


# These commands copy your files into the specified directory in the image
# and set that as the working location
COPY . /usr/src/gridsimplex
WORKDIR /usr/src/gridsimplex

# This command compiles your app using GCC, adjust for your source code
RUN mpiCC -o myapp main.cpp

# This command runs your application, comment out this line to compile only
CMD ["./myapp"]

LABEL Name=gridsimplex Version=0.0.1
