FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        autodock-vina \
        bash \
        ca-certificates \
        openbabel \
        pymol \
        python3 \
        xvfb \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work

CMD ["bash", "scripts/run_workflow.sh"]
