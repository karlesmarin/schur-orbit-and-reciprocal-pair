#!/bin/sh
# Copyright (c) 2026 Carles Marín. All rights reserved.
# Carles Marín <karlesmarin@gmail.com>, con Claude (Anthropic) como asistente.
# usage: sh sg.sh script.sage [args...]   -- runs sage in docker on the directory of this file
HERE=$(cd "$(dirname "$0")" && (pwd -W 2>/dev/null || pwd))
MSYS_NO_PATHCONV=1 docker run --rm -v "$HERE:/work" -w /work sage-normaliz:local sage "$@"
