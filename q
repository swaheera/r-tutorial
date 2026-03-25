library(dplyr)
library(forecast)

# ==============================================================================
# SCENARIOS & TIME PERIODS
# ==============================================================================
scenario_pcts <- c(-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3)

# ==============================================================================
# TIME WINDOWS — manually written
#
#   Training always starts April 2020.
#   Test is always 12 months.
#   Last data point is January 2026.
#
#   Window  1: train Apr 2020 – Mar 2024,  test Apr 2024 – Mar 2025
#   Window  2: train Apr 2020 – Apr 2024,  test May 2024 – Apr 2025
#   Window  3: train Apr 2020 – May 2024,  test Jun 2024 – May 2025
#   Window  4: train Apr 2020 – Jun 2024,  test Jul 2024 – Jun 2025
#   Window  5: train Apr 2020 – Jul 2024,  test Aug 2024 – Jul 2025
#   Window  6: train Apr 2020 – Aug 2024,  test Sep 2024 – Aug 2025
#   Window  7: train Apr 2020 – Sep 2024,  test Oct 2024 – Sep 2025
#   Window  8: train Apr 2020 – Oct 2024,  test Nov 2024 – Oct 2025
#   Window  9: train Apr 2020 – Nov 2024,  test Dec 2024 – Nov 2025
#   Window 10: train Apr 2020 – Dec 2024,  test Jan 2025 – Dec 2025
#   Window 11: train Apr 2020 – Jan 2025,  test Feb 2025 – Jan 2026
# ==============================================================================
windows <- list(
  list(train_start = "2020-04-01", train_end = "2024-03-31", test_start = "2024-04-01", test_end = "2025-03-31"),
  list(train_start = "2020-04-01", train_end = "2024-04-30", test_start = "2024-05-01", test_end = "2025-04-30"),
  list(train_start = "2020-04-01", train_end = "2024-05-31", test_start = "2024-06-01", test_end = "2025-05-31"),
  list(train_start = "2020-04-01", train_end = "2024-06-30", test_start = "2024-07-01", test_end = "2025-06-30"),
  list(train_start = "2020-04-01", train_end = "2024-07-31", test_start = "2024-08-01", test_end = "2025-07-31"),
  list(train_start = "2020-04-01", train_end = "2024-08-31", test_start = "2024-09-01", test_end = "2025-08-31"),
  list(train_start = "2020-04-01", train_end = "2024-09-30", test_start = "2024-10-01", test_end = "2025-09-30"),
  list(train_start = "2020-04-01", train_end = "2024-10-31", test_start = "2024-11-01", test_end = "2025-10-31"),
  list(train_start = "2020-04-01", train_end = "2024-11-30", test_start = "2024-12-01", test_end = "2025-11-30"),
  list(train_start = "2020-04-01", train_end = "2024-12-31", test_start = "2025-01-01", test_end = "2025-12-31"),
  list(train_start = "2020-04-01", train_end = "2025-01-31", test_start = "2025-02-01", test_end = "2026-01-31")
)

cat("Number of time windows:", length(windows), "\n")
cat("Number of scenarios:   ", length(scenario_pcts), "\n")
cat("Total runs:            ", length(windows) * length(scenario_pcts), "\n\n")

# ==============================================================================
# PREP — lag population once (doesn't change across runs)
# ==============================================================================
df <- mydf %>%
  arrange(date) %>%
  mutate(pop_lag1 = dplyr::lag(population, 1))

# ==============================================================================
# MAIN LOOP — simple nested for-loops, no fancy functions
# ==============================================================================
all_results <- data.frame()

for (w in seq_along(windows)) {

  win <- windows[[w]]

  train <- df %>%
    dplyr::filter(date >= win$train_start & date <= win$train_end & !is.na(pop_lag1))
  test  <- df %>%
    dplyr::filter(date >= win$test_start  & date <= win$test_end  & !is.na(pop_lag1))

  if (nrow(train) == 0 | nrow(test) == 0) next

  # Fit ARIMA models once per window (training data doesn't change across scenarios)
  xreg_train <- log(train$pop_lag1)

  model_total  <- auto.arima(log(train$total),  xreg = xreg_train)
  model_apples <- auto.arima(log(train$apples), xreg = xreg_train)
  model_banana <- auto.arima(log(train$banana), xreg = xreg_train)

  cat("Window", w, ":", win$train_start, "to", win$train_end,
      "| test:", win$test_start, "to", win$test_end, "\n")

  for (scenario_pct in scenario_pcts) {

    test_run <- test  # fresh copy

    # --- scenario population ---
    k        <- nrow(test_run)
    base_pop <- train$population[nrow(train)]
    end_pop  <- base_pop * (1 + scenario_pct / 100)

    test_run$pop_under_scenario     <- base_pop + (1:k) / k * (end_pop - base_pop)
    test_run$log_pop_under_scenario <- log(test_run$pop_under_scenario)

    # --- forecast ---
    xreg_test <- test_run$log_pop_under_scenario

    fc_total  <- forecast(model_total,  xreg = xreg_test)
    fc_apples <- forecast(model_apples, xreg = xreg_test)
    fc_banana <- forecast(model_banana, xreg = xreg_test)

    test_run$fc_total          <- exp(as.numeric(fc_total$mean))
    test_run$fc_apples         <- exp(as.numeric(fc_apples$mean))
    test_run$fc_banana_method1 <- exp(as.numeric(fc_banana$mean))
    test_run$fc_banana_method2 <- test_run$fc_total - test_run$fc_apples

    # --- log forecasts ---
    test_run$log_fc_total          <- as.numeric(fc_total$mean)
    test_run$log_fc_apples         <- as.numeric(fc_apples$mean)
    test_run$log_fc_banana_method1 <- as.numeric(fc_banana$mean)
    test_run$log_fc_banana_method2 <- log(test_run$fc_banana_method2)

    # --- cumulative absolute error (%) ---
    test_run$cum_err_total_pct   <- abs(cumsum(test_run$total)  - cumsum(test_run$fc_total))          / cumsum(test_run$total)  * 100
    test_run$cum_err_apples_pct  <- abs(cumsum(test_run$apples) - cumsum(test_run$fc_apples))         / cumsum(test_run$apples) * 100
    test_run$cum_err_banana1_pct <- abs(cumsum(test_run$banana) - cumsum(test_run$fc_banana_method1)) / cumsum(test_run$banana) * 100
    test_run$cum_err_banana2_pct <- abs(cumsum(test_run$banana) - cumsum(test_run$fc_banana_method2)) / cumsum(test_run$banana) * 100

    # --- true population growth ---
    test_run$true_pop_growth_rate_pct <- (test_run$population / test_run$population[1] - 1) * 100

    # --- helper values for formula strings ---
    cum_actual_total  <- cumsum(test_run$total)
    cum_fc_total      <- cumsum(test_run$fc_total)
    cum_actual_apples <- cumsum(test_run$apples)
    cum_fc_apples     <- cumsum(test_run$fc_apples)
    cum_actual_banana <- cumsum(test_run$banana)
    cum_fc_banana1    <- cumsum(test_run$fc_banana_method1)
    cum_fc_banana2    <- cumsum(test_run$fc_banana_method2)
    pop_first         <- test_run$population[1]

    # --- build results row(s) with % signs and formula columns ---
    run_results <- test_run %>%
      dplyr::transmute(
        training_period                = paste(win$train_start, "to", win$train_end),
        testing_period                 = paste(win$test_start,  "to", win$test_end),
        date                           = date,
        scenario                       = paste0(scenario_pct, "%"),

        base_pop_scenario              = base_pop,
        base_pop_date                  = train$date[nrow(train)],

        pop_under_scenario             = pop_under_scenario,
        pop_under_scenario_formula     = paste0(round(base_pop,2), " + ", row_number(), "/", k,
                                                " * (", round(end_pop,2), " - ", round(base_pop,2), ")"),

        log_pop_under_scenario         = log_pop_under_scenario,
        log_pop_under_scenario_formula = paste0("log(", round(pop_under_scenario, 2), ")"),

        true_pop                       = population,

        true_pop_growth_rate_pct       = paste0(round(true_pop_growth_rate_pct, 4), "%"),
        true_pop_growth_rate_formula   = paste0("(", population, " / ", pop_first, " - 1) * 100%"),

        true_total                     = total,
        true_apple                     = apples,
        true_banana                    = banana,

        total_forecast                 = fc_total,
        total_forecast_formula         = paste0("exp(", round(log_fc_total, 6), ")"),

        apple_forecast                 = fc_apples,
        apple_forecast_formula         = paste0("exp(", round(log_fc_apples, 6), ")"),

        banana_method1_forecast        = fc_banana_method1,
        banana_method1_forecast_formula = paste0("exp(", round(log_fc_banana_method1, 6), ")"),

        banana_method2_forecast        = fc_banana_method2,
        banana_method2_forecast_formula = paste0(round(fc_total, 4), " - ", round(fc_apples, 4)),

        log_total_forecast             = log_fc_total,
        log_total_forecast_formula     = paste0("log(", round(fc_total, 4), ")"),

        log_apple_forecast             = log_fc_apples,
        log_apple_forecast_formula     = paste0("log(", round(fc_apples, 4), ")"),

        log_banana_method1_forecast    = log_fc_banana_method1,
        log_banana_method1_forecast_formula = paste0("log(", round(fc_banana_method1, 4), ")"),

        log_banana_method2_forecast    = log_fc_banana_method2,
        log_banana_method2_forecast_formula = paste0("log(", round(fc_banana_method2, 4), ")"),

        cum_error_total_pct            = paste0(round(cum_err_total_pct, 4), "%"),
        cum_error_total_formula        = paste0("abs(", round(cum_actual_total, 2), " - ",
                                                round(cum_fc_total, 2), ") / ",
                                                round(cum_actual_total, 2), " * 100%"),

        cum_error_apple_pct            = paste0(round(cum_err_apples_pct, 4), "%"),
        cum_error_apple_formula        = paste0("abs(", round(cum_actual_apples, 2), " - ",
                                                round(cum_fc_apples, 2), ") / ",
                                                round(cum_actual_apples, 2), " * 100%"),

        cum_error_banana_method1_pct   = paste0(round(cum_err_banana1_pct, 4), "%"),
        cum_error_banana_method1_formula = paste0("abs(", round(cum_actual_banana, 2), " - ",
                                                  round(cum_fc_banana1, 2), ") / ",
                                                  round(cum_actual_banana, 2), " * 100%"),

        cum_error_banana_method2_pct   = paste0(round(cum_err_banana2_pct, 4), "%"),
        cum_error_banana_method2_formula = paste0("abs(", round(cum_actual_banana, 2), " - ",
                                                  round(cum_fc_banana2, 2), ") / ",
                                                  round(cum_actual_banana, 2), " * 100%")
      )

    all_results <- rbind(all_results, run_results)
  }
}

cat("\nDone! Total rows:", nrow(all_results), "\n")
print(head(all_results, 20))
