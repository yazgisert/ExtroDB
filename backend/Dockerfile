FROM python:3.12-slim

# Set the working directory inside the container
WORKDIR /app

# Install Python dependencies first (separate step so Docker can cache this
# layer — if only the app code changes, this layer doesn't need to rebuild)
COPY backend/requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy only the backend application code
COPY scripts/indexing/searchable_resource_webapp.py .
COPY searchable_resource/ searchable_resource/

# Expose the port uvicorn will listen on
EXPOSE 8000

# Start the server
CMD ["uvicorn", "searchable_resource_webapp:app", "--host", "0.0.0.0", "--port", "8000"]
