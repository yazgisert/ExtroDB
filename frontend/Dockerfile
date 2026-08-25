FROM nginx:alpine

RUN rm -rf /usr/share/nginx/html/*

COPY searchable_resource_app/ /usr/share/nginx/html/
COPY frontend/nginx.conf /etc/nginx/conf.d/default.conf

RUN chmod -R 755 /usr/share/nginx/html

EXPOSE 80
