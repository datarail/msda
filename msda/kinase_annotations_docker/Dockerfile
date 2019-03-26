# Use an official Python runtime as a parent image
FROM python:2.7

# Copy the current directory contents into the container at /app
COPY . /app
  
# Install any needed packages specified in requirements.txt
RUN cd /app \
    && pip install --trusted-host pypi.python.org -r requirements.txt \
    && chmod u+x app.py \
    && cc -O3 -o ./resources/NetPhorest/netphorest ./resources/NetPhorest/netphorest.c -lm


# Set the working directory to /app
WORKDIR /app

VOLUME /input
VOLUME /output

# Run app.py when the container launches
CMD ["python", "app.py"]