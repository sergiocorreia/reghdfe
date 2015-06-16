// -------------------------------------------------------------
// Simple assertions
// -------------------------------------------------------------
cap pr drop Assert
program define Assert
    syntax anything(everything equalok) [, MSG(string asis) RC(integer 198)]
    if !(`anything') {
        di as error `msg'
        exit `rc'
    }
end
