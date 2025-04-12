using TestItemRunner
using Aqua

@testitem "Look at method ambiguities" begin
    using Aqua
    # Disabling dependency compatibility because it's giving spurious output.
    Aqua.test_all(
        BijectiveHilbert;
        deps_compat=false
        )
end
