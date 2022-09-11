####################
# Input versioning #
####################

const INPUT_VERSIONS = [
    [v"0.15"],
]

function same_minor(v1::VersionNumber, v2::VersionNumber)
    Base.thisminor(v1) == Base.thisminor(v2)
end

function find_version_group(version::VersionNumber)
    findfirst(group -> any(v -> same_minor(v, version), group), INPUT_VERSIONS)
end

function compatible_input_versions(version::VersionNumber)
    i = find_version_group(version)
    i === nothing && return nothing
    map(INPUT_VERSIONS[i]) do v
        replace(string(v), r"\.0$" => "")
    end
end

function check_input_version(version::VersionNumber)
    i = find_version_group(version)
    j = find_version_group(PACKAGE_VERSION)
    j === nothing && error("Unreachable")
    if i === nothing
        error(
            "Input file v$version is not supported. " *
            "Available input file versions for MarbleRun v$PACKAGE_VERSION are $(compatible_input_versions(PACKAGE_VERSION))."
        )
    end
    if i != j
        error(
            "Input file v$version is not compatible with MarbleRun v$PACKAGE_VERSION. " *
            "Available MarbleRun versions are $(compatible_input_versions(version))."
        )
    end
end
