# Self-contained MinIO image for htslib's S3 integration tests.
#
# Bundles the MinIO server with the `mc` client and an entrypoint that
# provisions a fixed test bucket on startup, so the container is ready to
# receive S3 traffic (via hfile_s3.c) as soon as its healthcheck passes.
#
# Build:
#   docker build -t htslib-minio-test test/docker/minio
# Run:
#   docker run --rm -p 9000:9000 \
#     -e MINIO_ROOT_USER=minioadmin -e MINIO_ROOT_PASSWORD=minioadmin \
#     -e MINIO_TEST_BUCKET=htslib-test \
#     htslib-minio-test

FROM minio/mc:latest AS mc

FROM minio/minio:latest

COPY --from=mc /usr/bin/mc /usr/bin/mc
COPY bootstrap.sh /usr/local/bin/bootstrap.sh

ENV MINIO_ROOT_USER=minioadmin \
    MINIO_ROOT_PASSWORD=minioadmin \
    MINIO_TEST_BUCKET=htslib-test

EXPOSE 9000 9001

ENTRYPOINT ["/usr/local/bin/bootstrap.sh"]
