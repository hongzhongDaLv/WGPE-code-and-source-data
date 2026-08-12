# Configuration

Copy `config.example.yml` to `config.yml`. Paths may be absolute or relative to the repository root. `project_root` must be a dedicated working directory because a full rerun writes derived products beneath it. The preparation helper validates every required file, reconstructs the original script layout in that working directory and writes the machine-readable analysis configuration. Users do not edit individual analysis scripts.
