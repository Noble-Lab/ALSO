
# #  Take the AS max for each fragment ------------------------------------------
# #+ 
# #+ Taking the AS max (pmax()) gives us the best AS score per read pair
# 
# #  mm10
# command <- paste0(
#     "<- ", "pmax(",
#         variable_uniline, "$AS.mm10.odd, ",
#         variable_uniline, "$AS.mm10.even",
#     ")"
# )
# operation <- makeOperation(paste0(variable_uniline, "$AS.mm10.pmax"), command)
# evaluateOperation(operation)
# 
# #  129S1
# command <- paste0(
#     "<- ", "pmax(",
#         variable_uniline, "$AS.129S1.odd, ",
#         variable_uniline, "$AS.129S1.even",
#     ")"
# )
# operation <- makeOperation(paste0(variable_uniline, "$AS.129S1.pmax"), command)
# evaluateOperation(operation)
# 
# #  CAST
# command <- paste0(
#     "<- ", "pmax(",
#         variable_uniline, "$AS.CAST.odd, ",
#         variable_uniline, "$AS.CAST.even",
#     ")"
# )
# operation <- makeOperation(paste0(variable_uniline, "$AS.CAST.pmax"), command)
# evaluateOperation(operation)

# AS.pmax <- uniline.m.full %>%
#     dplyr::select(
#         coordinate.odd, coordinate.even, AS.mm10.pmax, AS.129S1.pmax, AS.CAST.pmax
#     )
# colnames(AS.pmax) <- gsub(".pmax", "", colnames(AS.pmax))

# AS.pmax <- AS.pmax %>% 
#     dplyr::mutate(difference = AS.129S1 - AS.CAST) %>% 
#     dplyr::mutate(
#         assignment.initial = case_when(
#             difference >= (-1 * int) & difference <= int ~ "Neutral",
#             difference > int ~ "129S1-SvImJ",
#             difference < (-1 * int) ~ "CAST-EiJ"
#         )
#     )

# AS.pmax$assignment <- ifelse(
#     is.na(AS.pmax$assignment.initial),
#     ifelse(
#         !is.na(AS.pmax$AS.129S1),
#         "129S1-SvImJ",
#         ifelse(
#             !is.na(AS.pmax$AS.CAST),
#             "CAST-EiJ",
#             AS.pmax$assignment.initial
#         )
#     ),
#     AS.pmax$assignment.initial
# ) %>%
#     forcats::as_factor()


# AS.pmax <- AS.pmax %>%
#     dplyr::select(
#         coordinate.odd, coordinate.even, AS.mm10, AS.129S1, AS.CAST, assignment
#     )

# AS.pmax$ref.mm10 <- ifelse(is.na(AS.pmax$AS.mm10), "0", "1")
# AS.pmax$ref.129S1 <- ifelse(is.na(AS.pmax$AS.129S1), "0", "1")
# AS.pmax$ref.CAST <- ifelse(is.na(AS.pmax$AS.CAST), "0", "1")
# 
# AS.pmax$trinary <- paste0(
#     AS.pmax$ref.mm10,
#     AS.pmax$ref.129S1,
#     AS.pmax$ref.CAST
# ) %>%
#     forcats::as_factor()
# 
# AS.pmax$assignment_trinary <- paste(
#     AS.pmax$assignment,
#     AS.pmax$trinary
# ) %>%
#     forcats::as_factor()
# 
# AS.pmax$trinary.r <- AS.pmax$trinary
# AS.pmax$assignment_trinary.r <- AS.pmax$assignment_trinary
# 
# AS.pmax$trinary.r <- AS.pmax$trinary.r %>% plyr::revalue(  # old_name = new_name
#     .,
#     c(
#         "100" = "N-masked mm10 reference",
#         "111" = "All three references",
#         "110" = "N-masked mm10 and 129 references",
#         "101" = "N-masked mm10 and CAST references",
#         "011" = "129 and CAST references",
#         "010" = "129 reference",
#         "001" = "CAST reference"
#     )
# )
# AS.pmax$assignment_trinary.r <- AS.pmax$assignment_trinary.r %>% plyr::revalue(
#     .,
#     c(
#         "NA 100" = "assignment: NA; alignment: N-masked mm10 reference only",
#         "Neutral 111" = "assignment: Neutral; alignment: all three references",
#         "129S1-SvImJ 111" = "assignment: 129; alignment: all three references",
#         "CAST-EiJ 111" = "assignment: CAST; alignment: all three references",
#         "129S1-SvImJ 110" = "assignment: 129; alignment: N-masked mm10 and 129 references",
#         "CAST-EiJ 101" = "assignment: CAST; alignment: N-masked mm10 and CAST references",
#         "CAST-EiJ 011" = "assignment: CAST; alignment: 129 and CAST references",
#         "129S1-SvImJ 010" = "assignment: 129; alignment: 129 reference only",
#         "Neutral 011" = "assignment: Neutral; alignment: 129 and CAST references",
#         "129S1-SvImJ 011" = "assignment: 129; alignment: 129 and CAST references",
#         "CAST-EiJ 001" = "assignment: CAST; alignment: CAST reference only"
#     )
# )

#  More stuff, TBD ------------------------------------------------------------
# AS.pmax <- AS.pmax %>%
#     select(
#         coordinate.odd,
#         coordinate.even,
#         AS.mm10,
#         AS.129S1,
#         AS.CAST,
#         assignment,
#         trinary,
#         trinary.r,
#         assignment_trinary,
#         assignment_trinary.r
#     )
