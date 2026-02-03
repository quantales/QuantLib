/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2010 Klaus Spanderen
 Copyright (C) 2025 Yassine Idyiahia
 
 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <https://www.quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include <ql/math/functional.hpp>
#include <ql/pricingengines/basket/kirkengine.hpp>
#include <ql/pricingengines/blackcalculator.hpp>
#include <utility>

namespace QuantLib {

    KirkEngine::KirkEngine(ext::shared_ptr<GeneralizedBlackScholesProcess> process1,
                           ext::shared_ptr<GeneralizedBlackScholesProcess> process2,
                           Real correlation,
                           bool useSkewCorrection)
    : SpreadBlackScholesVanillaEngine(std::move(process1), std::move(process2), correlation),
      useSkewCorrection_(useSkewCorrection) {
    }

    Real KirkEngine::calculate(
        Real f1, Real f2, Real strike, Option::Type optionType,
        Real variance1, Real variance2, DiscountFactor df) const {

        const Real K_eff = f2 + strike;
        const Real f = f1 / K_eff;
        const Real w = f2 / K_eff;

        const Real sigma1 = std::sqrt(variance1);
        const Real sigma2 = std::sqrt(variance2);

        // Kirk's volatility
        const Real sigma_kirk_sq = variance1
                                 + variance2 * squared(w)
                                 - 2.0 * rho_ * sigma1 * sigma2 * w;

        Real v = std::sqrt(std::max(sigma_kirk_sq, 0.0));

        if (useSkewCorrection_ && v > QL_EPSILON) {
            // Alos-Leon (2015) skew correction:
            // sigma_modified = sigma_kirk + skew_slope * (log(f1) - log(K_eff))
            const Real sigma_kirk = v;
            const Real numerator_term = squared(sigma2 * w - rho_ * sigma1);
            const Real skew_factor = variance2 * f2 * strike / squared(K_eff);
            const Real sigma_kirk_cubed = sigma_kirk * sigma_kirk * sigma_kirk;
            const Real skew_slope = 0.5 * numerator_term * skew_factor / sigma_kirk_cubed;

            const Real x = std::log(f1);
            const Real x_atm = std::log(K_eff);

            v = sigma_kirk + skew_slope * (x - x_atm);

            // Ensure corrected volatility remains positive
            if (v <= 0.0) {
                v = sigma_kirk;
            }
        }

        BlackCalculator black(
             ext::make_shared<PlainVanillaPayoff>(optionType, 1.0),
             f, v, df);

        return K_eff * black.value();
    }
}

