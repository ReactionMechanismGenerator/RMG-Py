tree(
"""
L1: R
    L2: C
        L3: Cbf
            L4: Cbf-CbCbCbf
            L4: Cbf-CbCbfCbf
            L4: Cbf-CbfCbfCbf
        L3: Cb
            L4: Cb-H
            L4: Cb-Os
            L4: Cb-C
                L5: Cb-Cs
                L5: Cb-Cds
                    L6: Cb-(Cds-Od)
                    L6: Cb-(Cds-Cd)
                        L7: Cb-(Cds-Cds)
                        L7: Cb-(Cds-Cdd)
                            L8: Cb-(Cds-Cdd-Od)
                            L8: Cb-(Cds-Cdd-Cd)
                L5: Cb-Ct
                L5: Cb-Cb
        L3: Ct
            L4: Ct-H
            L4: Ct-Os
            L4: Ct-C
                L5: Ct-Cs
                L5: Ct-Cds
                    L6: Ct-(Cds-Od)
                    L6: Ct-(Cds-Cd)
                        L7: Ct-(Cds-Cds)
                        L7: Ct-(Cds-Cdd)
                            L8: Ct-(Cds-Cdd-Od)
                            L8: Ct-(Cds-Cdd-Cd)
                L5: Ct-Ct
                L5: Ct-Cb
        L3: Cdd
            L4: Cdd-OdOd
            L4: Cdd-CdOd
                L5: Cdd-CdsOd
                L5: Cdd-CddOd
                    L6: Cdd-(Cdd-Od)Od
                    L6: Cdd-(Cdd-Cd)Od
            L4: Cdd-CdCd
                L5: Cdd-CddCdd
                    L6: Cdd-(Cdd-Od)(Cdd-Od)
                    L6: Cdd-(Cdd-Od)(Cdd-Cd)
                    L6: Cdd-(Cdd-Cd)(Cdd-Cd)
                L5: Cdd-CddCds
                    L6: Cdd-(Cdd-Od)Cds
                    L6: Cdd-(Cdd-Cd)Cds
                L5: Cdd-CdsCds
        L3: Cds
            L4: Cds-OdHH
            L4: Cds-OdOsH
            L4: Cds-OdOsOs
            L4: Cds-OdCH
                L5: Cds-OdCsH
                L5: Cds-OdCdsH
                    L6: Cds-Od(Cds-Od)H
                    L6: Cds-Od(Cds-Cd)H
                        L7: Cds-Od(Cds-Cds)H
                        L7: Cds-Od(Cds-Cdd)H
                            L8: Cds-Od(Cds-Cdd-Od)H
                            L8: Cds-Od(Cds-Cdd-Cd)H
                L5: Cds-OdCtH
                L5: Cds-OdCbH
            L4: Cds-OdCOs
                L5: Cds-OdCsOs
                L5: Cds-OdCdsOs
                    L6: Cds-Od(Cds-Od)Os
                    L6: Cds-Od(Cds-Cd)Os
                        L7: Cds-Od(Cds-Cds)Os
                        L7: Cds-Od(Cds-Cdd)Os
                            L8: Cds-Od(Cds-Cdd-Od)Os
                            L8: Cds-Od(Cds-Cdd-Cd)Os
                L5: Cds-OdCtOs
                L5: Cds-OdCbOs
            L4: Cds-OdCC
                L5: Cds-OdCsCs
                L5: Cds-OdCdsCs
                    L6: Cds-Od(Cds-Od)Cs
                    L6: Cds-Od(Cds-Cd)Cs
                        L7: Cds-Od(Cds-Cds)Cs
                        L7: Cds-Od(Cds-Cdd)Cs
                            L8: Cds-Od(Cds-Cdd-Od)Cs
                            L8: Cds-Od(Cds-Cdd-Cd)Cs
                L5: Cds-OdCdsCds
                    L6: Cds-Od(Cds-Od)(Cds-Od)
                    L6: Cds-Od(Cds-Cd)(Cds-Od)
                        L7: Cds-Od(Cds-Cds)(Cds-Od)
                        L7: Cds-Od(Cds-Cdd)(Cds-Od)
                            L8: Cds-Od(Cds-Cdd-Od)(Cds-Od)
                            L8: Cds-Od(Cds-Cdd-Cd)(Cds-Od)
                    L6: Cds-Od(Cds-Cd)(Cds-Cd)
                        L7: Cds-Od(Cds-Cds)(Cds-Cds)
                        L7: Cds-Od(Cds-Cdd)(Cds-Cds)
                            L8: Cds-Od(Cds-Cdd-Od)(Cds-Cds)
                            L8: Cds-Od(Cds-Cdd-Cd)(Cds-Cds)
                        L7: Cds-Od(Cds-Cdd)(Cds-Cdd)
                            L8: Cds-Od(Cds-Cdd-Od)(Cds-Cdd-Od)
                            L8: Cds-Od(Cds-Cdd-Cd)(Cds-Cdd-Od)
                            L8: Cds-Od(Cds-Cdd-Cd)(Cds-Cdd-Cd)
                L5: Cds-OdCtCs
                L5: Cds-OdCtCds
                    L6: Cds-OdCt(Cds-Od)
                    L6: Cds-OdCt(Cds-Cd)
                        L7: Cds-OdCt(Cds-Cds)
                        L7: Cds-OdCt(Cds-Cdd)
                            L8: Cds-OdCt(Cds-Cdd-Od)
                            L8: Cds-OdCt(Cds-Cdd-Cd)
                L5: Cds-OdCtCt
                L5: Cds-OdCbCs
                L5: Cds-OdCbCds
                    L6: Cds-OdCb(Cds-Od)
                    L6: Cds-OdCb(Cds-Cd)
                        L7: Cds-OdCb(Cds-Cds)
                        L7: Cds-OdCb(Cds-Cdd)
                            L8: Cds-OdCb(Cds-Cdd-Od)
                            L8: Cds-OdCb(Cds-Cdd-Cd)
                L5: Cds-OdCbCt
                L5: Cds-OdCbCb
            L4: Cds-CdHH
                L5: Cds-CdsHH
                L5: Cds-CddHH
                    L6: Cds-(Cdd-Od)HH
                    L6: Cds-(Cdd-Cd)HH
            L4: Cds-CdOsH
                L5: Cds-CdsOsH
                L5: Cds-CddOsH
                    L6: Cds-(Cdd-Od)OsH
                    L6: Cds-(Cdd-Cd)OsH
            L4: Cds-CdOsOs
                L5: Cds-CdsOsOs
                L5: Cds-CddOsOs
                    L6: Cds-(Cdd-Od)OsOs
                    L6: Cds-(Cdd-Cd)OsOs
            L4: Cds-CdCH
                L5: Cds-CdsCsH
                L5: Cds-CdsCdsH
                    L6: Cds-Cds(Cds-Od)H
                    L6: Cds-Cds(Cds-Cd)H
                        L7: Cds-Cds(Cds-Cds)H
                        L7: Cds-Cds(Cds-Cdd)H
                            L8: Cds-Cds(Cds-Cdd-Od)H
                            L8: Cds-Cds(Cds-Cdd-Cd)H
                L5: Cds-CdsCtH
                L5: Cds-CdsCbH
                L5: Cds-CddCsH
                    L6: Cds-(Cdd-Od)CsH
                    L6: Cds-(Cdd-Cd)CsH
                L5: Cds-CddCdsH
                    L6: Cds-(Cdd-Od)(Cds-Od)H
                    L6: Cds-(Cdd-Od)(Cds-Cd)H
                        L7: Cds-(Cdd-Od)(Cds-Cds)H
                        L7: Cds-(Cdd-Od)(Cds-Cdd)H
                            L8: Cds-(Cdd-Od)(Cds-Cdd-Od)H
                            L8: Cds-(Cdd-Od)(Cds-Cdd-Cd)H
                    L6: Cds-(Cdd-Cd)(Cds-Od)H
                    L6: Cds-(Cdd-Cd)(Cds-Cd)H
                        L7: Cds-(Cdd-Cd)(Cds-Cds)H
                        L7: Cds-(Cdd-Cd)(Cds-Cdd)H
                            L8: Cds-(Cdd-Cd)(Cds-Cdd-Od)H
                            L8: Cds-(Cdd-Cd)(Cds-Cdd-Cd)H
                L5: Cds-CddCtH
                    L6: Cds-(Cdd-Od)CtH
                    L6: Cds-(Cdd-Cd)CtH
                L5: Cds-CddCbH
                    L6: Cds-(Cdd-Od)CbH
                    L6: Cds-(Cdd-Cd)CbH
            L4: Cds-CdCO
                L5: Cds-CdsCsOs
                L5: Cds-CdsCdsOs
                    L6: Cds-Cds(Cds-Od)Os
                    L6: Cds-Cds(Cds-Cd)Os
                        L7: Cds-Cds(Cds-Cds)Os
                        L7: Cds-Cds(Cds-Cdd)Os
                            L8: Cds-Cds(Cds-Cdd-Od)Os
                            L8: Cds-Cds(Cds-Cdd-Cd)Os
                L5: Cds-CdsCtOs
                L5: Cds-CdsCbOs
                L5: Cds-CddCsOs
                    L6: Cds-(Cdd-Od)CsOs
                    L6: Cds-(Cdd-Cd)CsOs
                L5: Cds-CddCdsOs
                    L6: Cds-(Cdd-Od)(Cds-Od)Os
                    L6: Cds-(Cdd-Od)(Cds-Cd)Os
                        L7: Cds-(Cdd-Od)(Cds-Cds)Os
                        L7: Cds-(Cdd-Od)(Cds-Cdd)Os
                            L8: Cds-(Cdd-Od)(Cds-Cdd-Od)Os
                            L8: Cds-(Cdd-Od)(Cds-Cdd-Cd)Os
                    L6: Cds-(Cdd-Cd)(Cds-Cd)Os
                        L7: Cds-(Cdd-Cd)(Cds-Cds)Os
                        L7: Cds-(Cdd-Cd)(Cds-Cdd)Os
                            L8: Cds-(Cdd-Cd)(Cds-Cdd-Od)Os
                            L8: Cds-(Cdd-Cd)(Cds-Cdd-Cd)Os
                L5: Cds-CddCtOs
                    L6: Cds-(Cdd-Od)CtOs
                    L6: Cds-(Cdd-Cd)CtOs
                L5: Cds-CddCbOs
                    L6: Cds-(Cdd-Od)CbOs
                    L6: Cds-(Cdd-Cd)CbOs
            L4: Cds-CdCC
                L5: Cds-CdsCsCs
                L5: Cds-CdsCdsCs
                    L6: Cds-Cds(Cds-Od)Cs
                    L6: Cds-Cds(Cds-Cd)Cs
                        L7: Cds-Cds(Cds-Cds)Cs
                        L7: Cds-Cds(Cds-Cdd)Cs
                            L8: Cds-Cds(Cds-Cdd-Od)Cs
                            L8: Cds-Cds(Cds-Cdd-Cd)Cs
                L5: Cds-CdsCdsCds
                    L6: Cds-Cds(Cds-Od)(Cds-Od)
                    L6: Cds-Cds(Cds-Od)(Cds-Cd)
                        L7: Cds-Cds(Cds-Od)(Cds-Cds)
                        L7: Cds-Cds(Cds-Od)(Cds-Cdd)
                            L8: Cds-Cds(Cds-Od)(Cds-Cdd-Od)
                            L8: Cds-Cds(Cds-Od)(Cds-Cdd-Cd)
                    L6: Cds-Cds(Cds-Cd)(Cds-Cd)
                        L7: Cds-Cds(Cds-Cds)(Cds-Cds)
                        L7: Cds-Cds(Cds-Cds)(Cds-Cdd)
                            L8: Cds-Cds(Cds-Cds)(Cds-Cdd-Od)
                            L8: Cds-Cds(Cds-Cds)(Cds-Cdd-Cd)
                        L7: Cds-Cds(Cds-Cdd)(Cds-Cdd)
                            L8: Cds-Cds(Cds-Cdd-Od)(Cds-Cdd-Od)
                            L8: Cds-Cds(Cds-Cdd-Od)(Cds-Cdd-Cd)
                            L8: Cds-Cds(Cds-Cdd-Cd)(Cds-Cdd-Cd)
                L5: Cds-CdsCtCs
                L5: Cds-CdsCtCds
                    L6: Cds-CdsCt(Cds-Od)
                    L6: Cds-CdsCt(Cds-Cd)
                        L7: Cds-Cds(Cds-Cds)Ct
                        L7: Cds-Cds(Cds-Cdd)Ct
                            L8: Cds-Cds(Cds-Cdd-Od)Ct
                            L8: Cds-Cds(Cds-Cdd-Cd)Ct
                L5: Cds-CdsCtCt
                L5: Cds-CdsCbCs
                L5: Cds-CdsCbCds
                    L6: Cds-CdsCb(Cds-Od)
                    L6: Cds-Cds(Cds-Cd)Cb
                        L7: Cds-Cds(Cds-Cds)Cb
                        L7: Cds-Cds(Cds-Cdd)Cb
                            L8: Cds-Cds(Cds-Cdd-Od)Cb
                            L8: Cds-Cds(Cds-Cdd-Cd)Cb
                L5: Cds-CdsCbCt
                L5: Cds-CdsCbCb
                L5: Cds-CddCsCs
                    L6: Cds-(Cdd-Od)CsCs
                    L6: Cds-(Cdd-Cd)CsCs
                L5: Cds-CddCdsCs
                    L6: Cds-(Cdd-Od)(Cds-Od)Cs
                    L6: Cds-(Cdd-Od)(Cds-Cd)Cs
                        L7: Cds-(Cdd-Od)(Cds-Cds)Cs
                        L7: Cds-(Cdd-Od)(Cds-Cdd)Cs
                            L8: Cds-(Cdd-Od)(Cds-Cdd-Od)Cs
                            L8: Cds-(Cdd-Od)(Cds-Cdd-Cd)Cs
                    L6: Cds-(Cdd-Cd)(Cds-Cd)Cs
                        L7: Cds-(Cdd-Cd)(Cds-Cds)Cs
                        L7: Cds-(Cdd-Cd)(Cds-Cdd)Cs
                            L8: Cds-(Cdd-Cd)(Cds-Cdd-Od)Cs
                            L8: Cds-(Cdd-Cd)(Cds-Cdd-Cd)Cs
                L5: Cds-CddCdsCds
                    L6: Cds-(Cdd-Od)(Cds-Od)(Cds-Od)
                    L6: Cds-(Cdd-Od)(Cds-Cd)(Cds-Od)
                        L7: Cds-(Cdd-Od)(Cds-Cds)(Cds-Od)
                        L7: Cds-(Cdd-Od)(Cds-Cdd)(Cds-Od)
                            L8: Cds-(Cdd-Od)(Cds-Cdd-Od)(Cds-Od)
                            L8: Cds-(Cdd-Od)(Cds-Cdd-Cd)(Cds-Od)
                    L6: Cds-(Cdd-Od)(Cds-Cd)(Cds-Cd)
                        L7: Cds-(Cdd-Od)(Cds-Cds)(Cds-Cds)
                        L7: Cds-(Cdd-Od)(Cds-Cdd)(Cds-Cds)
                            L8: Cds-(Cdd-Od)(Cds-Cdd-Od)(Cds-Cds)
                            L8: Cds-(Cdd-Od)(Cds-Cdd-Cd)(Cds-Cds)
                        L7: Cds-(Cdd-Od)(Cds-Cdd)(Cds-Cdd)
                            L8: Cds-(Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)
                            L8: Cds-(Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)
                            L8: Cds-(Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)
                    L6: Cds-(Cdd-Cd)(Cds-Od)(Cds-Od)
                    L6: Cds-(Cdd-Cd)(Cds-Od)(Cds-Cd)
                        L7: Cds-(Cdd-Cd)(Cds-Od)(Cds-Cds)
                        L7: Cds-(Cdd-Cd)(Cds-Od)(Cds-Cdd)
                            L8: Cds-(Cdd-Cd)(Cds-Od)(Cds-Cdd-Od)
                            L8: Cds-(Cdd-Cd)(Cds-Od)(Cds-Cdd-Cd)
                    L6: Cds-(Cdd-Cd)(Cds-Cd)(Cds-Cd)
                        L7: Cds-(Cdd-Cd)(Cds-Cds)(Cds-Cds)
                        L7: Cds-(Cdd-Cd)(Cds-Cdd)(Cds-Cds)
                            L8: Cds-(Cdd-Cd)(Cds-Cdd-Od)(Cds-Cds)
                            L8: Cds-(Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cds)
                        L7: Cds-(Cdd-Cd)(Cds-Cdd)(Cds-Cdd)
                            L8: Cds-(Cdd-Cd)(Cds-Cdd-Od)(Cds-Cdd-Od)
                            L8: Cds-(Cdd-Cd)(Cds-Cdd-Od)(Cds-Cdd-Cd)
                            L8: Cds-(Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)
                L5: Cds-CddCtCs
                    L6: Cds-(Cdd-Od)CtCs
                    L6: Cds-(Cdd-Cd)CtCs
                L5: Cds-CddCtCds
                    L6: Cds-(Cdd-Od)(Cds-Od)Ct
                    L6: Cds-(Cdd-Od)(Cds-Cd)Ct
                        L7: Cds-(Cdd-Od)(Cds-Cds)Ct
                        L7: Cds-(Cdd-Od)(Cds-Cdd)Ct
                            L8: Cds-(Cdd-Od)(Cds-Cdd-Od)Ct
                            L8: Cds-(Cdd-Od)(Cds-Cdd-Cd)Ct
                    L6: Cds-(Cdd-Cd)(Cds-Cd)Ct
                        L7: Cds-(Cdd-Cd)(Cds-Cds)Ct
                        L7: Cds-(Cdd-Cd)(Cds-Cdd)Ct
                            L8: Cds-(Cdd-Cd)(Cds-Cdd-Od)Ct
                            L8: Cds-(Cdd-Cd)(Cds-Cdd-Cd)Ct
                L5: Cds-CddCtCt
                    L6: Cds-(Cdd-Od)CtCt
                    L6: Cds-(Cdd-Cd)CtCt
                L5: Cds-CddCbCs
                    L6: Cds-(Cdd-Od)CbCs
                    L6: Cds-(Cdd-Cd)CbCs
                L5: Cds-CddCbCds
                    L6: Cds-(Cdd-Od)(Cds-Od)Cb
                    L6: Cds-(Cdd-Od)(Cds-Cd)Cb
                        L7: Cds-(Cdd-Od)(Cds-Cds)Cb
                        L7: Cds-(Cdd-Od)(Cds-Cdd)Cb
                            L8: Cds-(Cdd-Od)(Cds-Cdd-Od)Cb
                            L8: Cds-(Cdd-Od)(Cds-Cdd-Cd)Cb
                    L6: Cds-(Cdd-Cd)(Cds-Cd)Cb
                        L7: Cds-(Cdd-Cd)(Cds-Cds)Cb
                        L7: Cds-(Cdd-Cd)(Cds-Cdd)Cb
                            L8: Cds-(Cdd-Cd)(Cds-Cdd-Od)Cb
                            L8: Cds-(Cdd-Cd)(Cds-Cdd-Cd)Cb
                L5: Cds-CddCbCt
                    L6: Cds-(Cdd-Od)CbCt
                    L6: Cds-(Cdd-Cd)CbCt
                L5: Cds-CddCbCb
                    L6: Cds-(Cdd-Od)CbCb
                    L6: Cds-(Cdd-Cd)CbCb
        L3: Cs
            L4: Cs-HHHH
            L4: Cs-CHHH
                L5: Cs-CsHHH
                L5: Cs-CdsHHH
                    L6: Cs-(Cds-Od)HHH
                    L6: Cs-(Cds-Cd)HHH
                        L7: Cs-(Cds-Cds)HHH
                        L7: Cs-(Cds-Cdd)HHH
                            L8: Cs-(Cds-Cdd-Od)HHH
                            L8: Cs-(Cds-Cdd-Cd)HHH
                L5: Cs-CtHHH
                L5: Cs-CbHHH
            L4: Cs-OsHHH
            L4: Cs-OsOsHH
            L4: Cs-OsOsOsH
            L4: Cs-CCHH
                L5: Cs-CsCsHH
                L5: Cs-CdsCsHH
                    L6: Cs-(Cds-Od)CsHH
                    L6: Cs-(Cds-Cd)CsHH
                        L7: Cs-(Cds-Cds)CsHH
                        L7: Cs-(Cds-Cdd)CsHH
                            L8: Cs-(Cds-Cdd-Od)CsHH
                            L8: Cs-(Cds-Cdd-Cd)CsHH
                L5: Cs-CdsCdsHH
                    L6: Cs-(Cds-Od)(Cds-Od)HH
                    L6: Cs-(Cds-Od)(Cds-Cd)HH
                        L7: Cs-(Cds-Od)(Cds-Cds)HH
                        L7: Cs-(Cds-Od)(Cds-Cdd)HH
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)HH
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)HH
                    L6: Cs-(Cds-Cd)(Cds-Cd)HH
                        L7: Cs-(Cds-Cds)(Cds-Cds)HH
                        L7: Cs-(Cds-Cdd)(Cds-Cds)HH
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cds)HH
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cds)HH
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)HH
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)HH
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)HH
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)HH
                L5: Cs-CtCsHH
                L5: Cs-CtCdsHH
                    L6: Cs-(Cds-Od)CtHH
                    L6: Cs-(Cds-Cd)CtHH
                        L7: Cs-(Cds-Cds)CtHH
                        L7: Cs-(Cds-Cdd)CtHH
                            L8: Cs-(Cds-Cdd-Od)CtHH
                            L8: Cs-(Cds-Cdd-Cd)CtHH
                L5: Cs-CtCtHH
                L5: Cs-CbCsHH
                L5: Cs-CbCdsHH
                    L6: Cs-(Cds-Od)CbHH
                    L6: Cs-(Cds-Cd)CbHH
                        L7: Cs-(Cds-Cds)CbHH
                        L7: Cs-(Cds-Cdd)CbHH
                            L8: Cs-(Cds-Cdd-Od)CbHH
                            L8: Cs-(Cds-Cdd-Cd)CbHH
                L5: Cs-CbCtHH
                L5: Cs-CbCbHH
            L4: Cs-CCCH
                L5: Cs-CsCsCsH
                L5: Cs-CdsCsCsH
                    L6: Cs-(Cds-Od)CsCsH
                    L6: Cs-(Cds-Cd)CsCsH
                        L7: Cs-(Cds-Cds)CsCsH
                        L7: Cs-(Cds-Cdd)CsCsH
                            L8: Cs-(Cds-Cdd-Od)CsCsH
                            L8: Cs-(Cds-Cdd-Cd)CsCsH
                L5: Cs-CtCsCsH
                L5: Cs-CbCsCsH
                L5: Cs-CdsCdsCsH
                    L6: Cs-(Cds-Od)(Cds-Od)CsH
                    L6: Cs-(Cds-Od)(Cds-Cd)CsH
                        L7: Cs-(Cds-Od)(Cds-Cds)CsH
                        L7: Cs-(Cds-Od)(Cds-Cdd)CsH
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)CsH
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)CsH
                    L6: Cs-(Cds-Cd)(Cds-Cd)CsH
                        L7: Cs-(Cds-Cds)(Cds-Cds)CsH
                        L7: Cs-(Cds-Cdd)(Cds-Cds)CsH
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cds)CsH
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cds)CsH
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)CsH
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CsH
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CsH
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CsH
                L5: Cs-CtCdsCsH
                    L6: Cs-(Cds-Od)CtCsH
                    L6: Cs-(Cds-Cd)CtCsH
                        L7: Cs-(Cds-Cds)CtCsH
                        L7: Cs-(Cds-Cdd)CtCsH
                            L8: Cs-(Cds-Cdd-Od)CtCsH
                            L8: Cs-(Cds-Cdd-Cd)CtCsH
                L5: Cs-CbCdsCsH
                    L6: Cs-(Cds-Od)CbCsH
                    L6: Cs-(Cds-Cd)CbCsH
                        L7: Cs-(Cds-Cds)CbCsH
                        L7: Cs-(Cds-Cdd)CbCsH
                            L8: Cs-(Cds-Cdd-Od)CbCsH
                            L8: Cs-(Cds-Cdd-Cd)CbCsH
                L5: Cs-CtCtCsH
                L5: Cs-CbCtCsH
                L5: Cs-CbCbCsH
                L5: Cs-CdsCdsCdsH
                    L6: Cs-(Cds-Od)(Cds-Od)(Cds-Od)H
                    L6: Cs-(Cds-Od)(Cds-Od)(Cds-Cd)H
                        L7: Cs-(Cds-Od)(Cds-Od)(Cds-Cds)H
                        L7: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd)H
                            L8: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)H
                            L8: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)H
                    L6: Cs-(Cds-Od)(Cds-Cd)(Cds-Cd)H
                        L7: Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)H
                        L7: Cs-(Cds-Od)(Cds-Cdd)(Cds-Cds)H
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)H
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cds)H
                        L7: Cs-(Cds-Od)(Cds-Cdd)(Cds-Cdd)H
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)H
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)H
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)H
                    L6: Cs-(Cds-Cd)(Cds-Cd)(Cds-Cd)H
                        L7: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H
                        L7: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)H
                            L8: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)H
                            L8: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)H
                        L7: Cs-(Cds-Cds)(Cds-Cdd)(Cds-Cdd)H
                            L8: Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)H
                            L8: Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Cd)H
                            L8: Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)H
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)H
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)H
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)H
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)H
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)H
                L5: Cs-CtCdsCdsH
                    L6: Cs-(Cds-Od)(Cds-Od)CtH
                    L6: Cs-(Cds-Od)(Cds-Cd)CtH
                        L7: Cs-(Cds-Od)(Cds-Cds)CtH
                        L7: Cs-(Cds-Od)(Cds-Cdd)CtH
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)CtH
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)CtH
                    L6: Cs-(Cds-Cd)(Cds-Cd)CtH
                        L7: Cs-(Cds-Cds)(Cds-Cds)CtH
                        L7: Cs-(Cds-Cdd)(Cds-Cds)CtH
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cds)CtH
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cds)CtH
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)CtH
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CtH
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CtH
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtH
                L5: Cs-CbCdsCdsH
                    L6: Cs-(Cds-Od)(Cds-Od)CbH
                    L6: Cs-(Cds-Od)(Cds-Cd)CbH
                        L7: Cs-(Cds-Od)(Cds-Cds)CbH
                        L7: Cs-(Cds-Od)(Cds-Cdd)CbH
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)CbH
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)CbH
                    L6: Cs-(Cds-Cd)(Cds-Cd)CbH
                        L7: Cs-(Cds-Cds)(Cds-Cds)CbH
                        L7: Cs-(Cds-Cdd)(Cds-Cds)CbH
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cds)CbH
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cds)CbH
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)CbH
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CbH
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CbH
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbH
                L5: Cs-CtCtCdsH
                    L6: Cs-CtCt(Cds-Od)H
                    L6: Cs-CtCt(Cds-Cd)H
                        L7: Cs-CtCt(Cds-Cds)H
                        L7: Cs-CtCt(Cds-Cdd)H
                            L8: Cs-CtCt(Cds-Cdd-Od)H
                            L8: Cs-CtCt(Cds-Cdd-Cd)H
                L5: Cs-CbCtCdsH
                    L6: Cs-CbCt(Cds-Od)H
                    L6: Cs-CbCt(Cds-Cd)H
                        L7: Cs-CbCt(Cds-Cds)H
                        L7: Cs-CbCt(Cds-Cdd)H
                            L8: Cs-CbCt(Cds-Cdd-Od)H
                            L8: Cs-CbCt(Cds-Cdd-Cd)H
                L5: Cs-CbCbCdsH
                    L6: Cs-CbCb(Cds-Od)H
                    L6: Cs-CbCb(Cds-Cd)H
                        L7: Cs-CbCb(Cds-Cds)H
                        L7: Cs-CbCb(Cds-Cdd)H
                            L8: Cs-CbCb(Cds-Cdd-Od)H
                            L8: Cs-CbCb(Cds-Cdd-Cd)H
                L5: Cs-CtCtCtH
                L5: Cs-CbCtCtH
                L5: Cs-CbCbCtH
                L5: Cs-CbCbCbH
            L4: Cs-CCCC
                L5: Cs-CsCsCsCs
                L5: Cs-CdsCsCsCs
                    L6: Cs-(Cds-Od)CsCsCs
                    L6: Cs-(Cds-Cd)CsCsCs
                        L7: Cs-(Cds-Cds)CsCsCs
                        L7: Cs-(Cds-Cdd)CsCsCs
                            L8: Cs-(Cds-Cdd-Od)CsCsCs
                            L8: Cs-(Cds-Cdd-Cd)CsCsCs
                L5: Cs-CtCsCsCs
                L5: Cs-CbCsCsCs
                L5: Cs-CdsCdsCsCs
                    L6: Cs-(Cds-Od)(Cds-Od)CsCs
                    L6: Cs-(Cds-Od)(Cds-Cd)CsCs
                        L7: Cs-(Cds-Od)(Cds-Cds)CsCs
                        L7: Cs-(Cds-Od)(Cds-Cdd)CsCs
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)CsCs
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)CsCs
                    L6: Cs-(Cds-Cd)(Cds-Cd)CsCs
                        L7: Cs-(Cds-Cds)(Cds-Cds)CsCs
                        L7: Cs-(Cds-Cdd)(Cds-Cds)CsCs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cds)CsCs
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cds)CsCs
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)CsCs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CsCs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CsCs
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CsCs
                L5: Cs-CtCdsCsCs
                    L6: Cs-(Cds-Od)CtCsCs
                    L6: Cs-(Cds-Cd)CtCsCs
                        L7: Cs-(Cds-Cds)CtCsCs
                        L7: Cs-(Cds-Cdd)CtCsCs
                            L8: Cs-(Cds-Cdd-Od)CtCsCs
                            L8: Cs-(Cds-Cdd-Cd)CtCsCs
                L5: Cs-CbCdsCsCs
                    L6: Cs-(Cds-Od)CbCsCs
                    L6: Cs-(Cds-Cd)CbCsCs
                        L7: Cs-(Cds-Cds)CbCsCs
                        L7: Cs-(Cds-Cdd)CbCsCs
                            L8: Cs-(Cds-Cdd-Od)CbCsCs
                            L8: Cs-(Cds-Cdd-Cd)CbCsCs
                L5: Cs-CtCtCsCs
                L5: Cs-CbCtCsCs
                L5: Cs-CbCbCsCs
                L5: Cs-CdsCdsCdsCs
                    L6: Cs-(Cds-Od)(Cds-Od)(Cds-Od)Cs
                    L6: Cs-(Cds-Od)(Cds-Od)(Cds-Cd)Cs
                        L7: Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Cs
                        L7: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd)Cs
                            L8: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)Cs
                            L8: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)Cs
                    L6: Cs-(Cds-Od)(Cds-Cd)(Cds-Cd)Cs
                        L7: Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Cs
                        L7: Cs-(Cds-Od)(Cds-Cdd)(Cds-Cds)Cs
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)Cs
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cds)Cs
                        L7: Cs-(Cds-Od)(Cds-Cdd)(Cds-Cdd)Cs
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)Cs
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)Cs
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs
                    L6: Cs-(Cds-Cd)(Cds-Cd)(Cds-Cd)Cs
                        L7: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cs
                        L7: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)Cs
                            L8: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Cs
                            L8: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)Cs
                        L7: Cs-(Cds-Cds)(Cds-Cdd)(Cds-Cdd)Cs
                            L8: Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)Cs
                            L8: Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Cd)Cs
                            L8: Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)Cs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)Cs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)Cs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs
                L5: Cs-CtCdsCdsCs
                    L6: Cs-(Cds-Od)(Cds-Od)CtCs
                    L6: Cs-(Cds-Od)(Cds-Cd)CtCs
                        L7: Cs-(Cds-Od)(Cds-Cds)CtCs
                        L7: Cs-(Cds-Od)(Cds-Cdd)CtCs
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)CtCs
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)CtCs
                    L6: Cs-(Cds-Cd)(Cds-Cd)CtCs
                        L7: Cs-(Cds-Cds)(Cds-Cds)CtCs
                        L7: Cs-(Cds-Cdd)(Cds-Cds)CtCs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cds)CtCs
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cds)CtCs
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)CtCs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CtCs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CtCs
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtCs
                L5: Cs-CbCdsCdsCs
                    L6: Cs-(Cds-Od)(Cds-Od)CbCs
                    L6: Cs-(Cds-Od)(Cds-Cd)CbCs
                        L7: Cs-(Cds-Od)(Cds-Cds)CbCs
                        L7: Cs-(Cds-Od)(Cds-Cdd)CbCs
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)CbCs
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)CbCs
                    L6: Cs-(Cds-Cd)(Cds-Cd)CbCs
                        L7: Cs-(Cds-Cds)(Cds-Cds)CbCs
                        L7: Cs-(Cds-Cdd)(Cds-Cds)CbCs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cds)CbCs
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cds)CbCs
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)CbCs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CbCs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CbCs
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbCs
                L5: Cs-CtCtCdsCs
                    L6: Cs-(Cds-Od)CtCtCs
                    L6: Cs-(Cds-Cd)CtCtCs
                        L7: Cs-(Cds-Cds)CtCtCs
                        L7: Cs-(Cds-Cdd)CtCtCs
                            L8: Cs-(Cds-Cdd-Od)CtCtCs
                            L8: Cs-(Cds-Cdd-Cd)CtCtCs
                L5: Cs-CbCtCdsCs
                    L6: Cs-(Cds-Od)CbCtCs
                    L6: Cs-(Cds-Cd)CbCtCs
                        L7: Cs-(Cds-Cds)CbCtCs
                        L7: Cs-(Cds-Cdd)CbCtCs
                            L8: Cs-(Cds-Cdd-Od)CbCtCs
                            L8: Cs-(Cds-Cdd-Cd)CbCtCs
                L5: Cs-CbCbCdsCs
                    L6: Cs-(Cds-Od)CbCbCs
                    L6: Cs-(Cds-Cd)CbCbCs
                        L7: Cs-(Cds-Cds)CbCbCs
                        L7: Cs-(Cds-Cdd)CbCbCs
                            L8: Cs-(Cds-Cdd-Od)CbCbCs
                            L8: Cs-(Cds-Cdd-Cd)CbCbCs
                L5: Cs-CtCtCtCs
                L5: Cs-CbCtCtCs
                L5: Cs-CbCbCtCs
                L5: Cs-CbCbCbCs
                L5: Cs-CdsCdsCdsCds
                    L6: Cs-(Cds-Od)(Cds-Od)(Cds-Od)(Cds-Od)
                    L6: Cs-(Cds-Od)(Cds-Od)(Cds-Od)(Cds-Cd)
                        L7: Cs-(Cds-Od)(Cds-Od)(Cds-Od)(Cds-Cds)
                        L7: Cs-(Cds-Od)(Cds-Od)(Cds-Od)(Cds-Cdd)
                            L8: Cs-(Cds-Od)(Cds-Od)(Cds-Od)(Cds-Cdd-Od)
                            L8: Cs-(Cds-Od)(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)
                    L6: Cs-(Cds-Od)(Cds-Od)(Cds-Cd)(Cds-Cd)
                        L7: Cs-(Cds-Od)(Cds-Od)(Cds-Cds)(Cds-Cds)
                        L7: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd)(Cds-Cds)
                            L8: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)
                            L8: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)(Cds-Cds)
                        L7: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd)(Cds-Cdd)
                            L8: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)
                            L8: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)
                            L8: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)
                    L6: Cs-(Cds-Od)(Cds-Cd)(Cds-Cd)(Cds-Cd)
                        L7: Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cds)
                        L7: Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cdd)
                            L8: Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)
                            L8: Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)
                        L7: Cs-(Cds-Od)(Cds-Cds)(Cds-Cdd)(Cds-Cdd)
                            L8: Cs-(Cds-Od)(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)
                            L8: Cs-(Cds-Od)(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Cd)
                            L8: Cs-(Cds-Od)(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)
                        L7: Cs-(Cds-Od)(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)
                    L6: Cs-(Cds-Cd)(Cds-Cd)(Cds-Cd)(Cds-Cd)
                        L7: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)
                        L7: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd)
                            L8: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)
                            L8: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)
                        L7: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)(Cds-Cdd)
                            L8: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)
                            L8: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Cd)
                            L8: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)
                        L7: Cs-(Cds-Cds)(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)
                            L8: Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)
                            L8: Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)
                            L8: Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)
                            L8: Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)
                L5: Cs-CtCdsCdsCds
                    L6: Cs-(Cds-Od)(Cds-Od)(Cds-Od)Ct
                    L6: Cs-(Cds-Od)(Cds-Od)(Cds-Cd)Ct
                        L7: Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Ct
                        L7: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd)Ct
                            L8: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)Ct
                            L8: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)Ct
                    L6: Cs-(Cds-Od)(Cds-Cd)(Cds-Cd)Ct
                        L7: Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Ct
                        L7: Cs-(Cds-Od)(Cds-Cdd)(Cds-Cds)Ct
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)Ct
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cds)Ct
                        L7: Cs-(Cds-Od)(Cds-Cdd)(Cds-Cdd)Ct
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)Ct
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)Ct
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct
                    L6: Cs-(Cds-Cd)(Cds-Cd)(Cds-Cd)Ct
                        L7: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Ct
                        L7: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)Ct
                            L8: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Ct
                            L8: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)Ct
                        L7: Cs-(Cds-Cds)(Cds-Cdd)(Cds-Cdd)Ct
                            L8: Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)Ct
                            L8: Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Cd)Ct
                            L8: Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)Ct
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)Ct
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)Ct
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct
                L5: Cs-CbCdsCdsCds
                    L6: Cs-(Cds-Od)(Cds-Od)(Cds-Od)Cb
                    L6: Cs-(Cds-Od)(Cds-Od)(Cds-Cd)Cb
                        L7: Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Cb
                        L7: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd)Cb
                            L8: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)Cb
                            L8: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)Cb
                    L6: Cs-(Cds-Od)(Cds-Cd)(Cds-Cd)Cb
                        L7: Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Cb
                        L7: Cs-(Cds-Od)(Cds-Cdd)(Cds-Cds)Cb
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)Cb
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cds)Cb
                        L7: Cs-(Cds-Od)(Cds-Cdd)(Cds-Cdd)Cb
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)Cb
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)Cb
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb
                    L6: Cs-(Cds-Cd)(Cds-Cd)(Cds-Cd)Cb
                        L7: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cb
                        L7: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)Cb
                            L8: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Cb
                            L8: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)Cb
                        L7: Cs-(Cds-Cds)(Cds-Cdd)(Cds-Cdd)Cb
                            L8: Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)Cb
                            L8: Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Cd)Cb
                            L8: Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)Cb
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)Cb
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)Cb
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb
                L5: Cs-CtCtCdsCds
                    L6: Cs-(Cds-Od)(Cds-Od)CtCt
                    L6: Cs-(Cds-Od)(Cds-Cd)CtCt
                        L7: Cs-(Cds-Od)(Cds-Cds)CtCt
                        L7: Cs-(Cds-Od)(Cds-Cdd)CtCt
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)CtCt
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)CtCt
                    L6: Cs-(Cds-Cd)(Cds-Cd)CtCt
                        L7: Cs-(Cds-Cds)(Cds-Cds)CtCt
                        L7: Cs-(Cds-Cdd)(Cds-Cds)CtCt
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cds)CtCt
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cds)CtCt
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)CtCt
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CtCt
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CtCt
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtCt
                L5: Cs-CbCtCdsCds
                    L6: Cs-(Cds-Od)(Cds-Od)CbCt
                    L6: Cs-(Cds-Od)(Cds-Cd)CbCt
                        L7: Cs-(Cds-Od)(Cds-Cds)CbCt
                        L7: Cs-(Cds-Od)(Cds-Cdd)CbCt
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)CbCt
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)CbCt
                    L6: Cs-(Cds-Cd)(Cds-Cd)CbCt
                        L7: Cs-(Cds-Cds)(Cds-Cds)CbCt
                        L7: Cs-(Cds-Cdd)(Cds-Cds)CbCt
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cds)CbCt
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cds)CbCt
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)CbCt
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CbCt
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CbCt
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbCt
                L5: Cs-CbCbCdsCds
                    L6: Cs-(Cds-Od)(Cds-Od)CbCb
                    L6: Cs-(Cds-Od)(Cds-Cd)CbCb
                        L7: Cs-(Cds-Od)(Cds-Cds)CbCb
                        L7: Cs-(Cds-Od)(Cds-Cdd)CbCb
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)CbCb
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)CbCb
                    L6: Cs-(Cds-Cd)(Cds-Cd)CbCb
                        L7: Cs-(Cds-Cds)(Cds-Cds)CbCb
                        L7: Cs-(Cds-Cdd)(Cds-Cds)CbCb
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cds)CbCb
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cds)CbCb
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)CbCb
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CbCb
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CbCb
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbCb
                L5: Cs-CtCtCtCds
                    L6: Cs-(Cds-Od)CtCtCt
                    L6: Cs-(Cds-Cd)CtCtCt
                        L7: Cs-(Cds-Cds)CtCtCt
                        L7: Cs-(Cds-Cdd)CtCtCt
                            L8: Cs-(Cds-Cdd-Od)CtCtCt
                            L8: Cs-(Cds-Cdd-Cd)CtCtCt
                L5: Cs-CbCtCtCds
                    L6: Cs-(Cds-Od)CbCtCt
                    L6: Cs-(Cds-Cd)CbCtCt
                        L7: Cs-(Cds-Cds)CbCtCt
                        L7: Cs-(Cds-Cdd)CbCtCt
                            L8: Cs-(Cds-Cdd-Od)CbCtCt
                            L8: Cs-(Cds-Cdd-Cd)CbCtCt
                L5: Cs-CbCbCtCds
                    L6: Cs-(Cds-Od)CbCbCt
                    L6: Cs-(Cds-Cd)CbCbCt
                        L7: Cs-(Cds-Cds)CbCbCt
                        L7: Cs-(Cds-Cdd)CbCbCt
                            L8: Cs-(Cds-Cdd-Od)CbCbCt
                            L8: Cs-(Cds-Cdd-Cd)CbCbCt
                L5: Cs-CbCbCbCds
                    L6: Cs-(Cds-Od)CbCbCb
                    L6: Cs-(Cds-Cd)CbCbCb
                        L7: Cs-(Cds-Cds)CbCbCb
                        L7: Cs-(Cds-Cdd)CbCbCb
                            L8: Cs-(Cds-Cdd-Od)CbCbCb
                            L8: Cs-(Cds-Cdd-Cd)CbCbCb
                L5: Cs-CtCtCtCt
                L5: Cs-CbCtCtCt
                L5: Cs-CbCbCtCt
                L5: Cs-CbCbCbCt
                L5: Cs-CbCbCbCb
            L4: Cs-CCCOs
                L5: Cs-CsCsCsOs
                L5: Cs-CdsCsCsOs
                    L6: Cs-(Cds-Od)CsCsOs
                    L6: Cs-(Cds-Cd)CsCsOs
                        L7: Cs-(Cds-Cds)CsCsOs
                        L7: Cs-(Cds-Cdd)CsCsOs
                            L8: Cs-(Cds-Cdd-Od)CsCsOs
                            L8: Cs-(Cds-Cdd-Cd)CsCsOs
                L5: Cs-OsCtCsCs
                L5: Cs-CbCsCsOs
                L5: Cs-CdsCdsCsOs
                    L6: Cs-(Cds-Od)(Cds-Od)CsOs
                    L6: Cs-(Cds-Od)(Cds-Cd)CsOs
                        L7: Cs-(Cds-Od)(Cds-Cds)CsOs
                        L7: Cs-(Cds-Od)(Cds-Cdd)CsOs
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)CsOs
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)CsOs
                    L6: Cs-(Cds-Cd)(Cds-Cd)CsOs
                        L7: Cs-(Cds-Cds)(Cds-Cds)CsOs
                        L7: Cs-(Cds-Cdd)(Cds-Cds)CsOs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cds)CsOs
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cds)CsOs
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)CsOs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CsOs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CsOs
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CsOs
                L5: Cs-CtCdsCsOs
                    L6: Cs-(Cds-Od)CtCsOs
                    L6: Cs-(Cds-Cd)CtCsOs
                        L7: Cs-(Cds-Cds)CtCsOs
                        L7: Cs-(Cds-Cdd)CtCsOs
                            L8: Cs-(Cds-Cdd-Od)CtCsOs
                            L8: Cs-(Cds-Cdd-Cd)CtCsOs
                L5: Cs-CbCdsCsOs
                    L6: Cs-(Cds-Od)CbCsOs
                    L6: Cs-(Cds-Cd)CbCsOs
                        L7: Cs-(Cds-Cds)CbCsOs
                        L7: Cs-(Cds-Cdd)CbCsOs
                            L8: Cs-(Cds-Cdd-Od)CbCsOs
                            L8: Cs-(Cds-Cdd-Cd)CbCsOs
                L5: Cs-CtCtCsOs
                L5: Cs-CbCtCsOs
                L5: Cs-CbCbCsOs
                L5: Cs-CdsCdsCdsOs
                    L6: Cs-(Cds-Od)(Cds-Od)(Cds-Od)Os
                    L6: Cs-(Cds-Od)(Cds-Od)(Cds-Cd)Os
                        L7: Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Os
                        L7: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd)Os
                            L8: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)Os
                            L8: Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)Os
                    L6: Cs-(Cds-Od)(Cds-Cd)(Cds-Cd)Os
                        L7: Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Os
                        L7: Cs-(Cds-Od)(Cds-Cdd)(Cds-Cds)Os
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)Os
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cds)Os
                        L7: Cs-(Cds-Od)(Cds-Cdd)(Cds-Cdd)Os
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)Os
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)Os
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Os
                    L6: Cs-(Cds-Cd)(Cds-Cd)(Cds-Cd)Os
                        L7: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Os
                        L7: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)Os
                            L8: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Os
                            L8: Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)Os
                        L7: Cs-(Cds-Cds)(Cds-Cdd)(Cds-Cdd)Os
                            L8: Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)Os
                            L8: Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Cd)Os
                            L8: Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Os
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)Os
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)Os
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)Os
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Os
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Os
                L5: Cs-CtCdsCdsOs
                    L6: Cs-(Cds-Od)(Cds-Od)CtOs
                    L6: Cs-(Cds-Od)(Cds-Cd)CtOs
                        L7: Cs-(Cds-Od)(Cds-Cds)CtOs
                        L7: Cs-(Cds-Od)(Cds-Cdd)CtOs
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)CtOs
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)CtOs
                    L6: Cs-(Cds-Cd)(Cds-Cd)CtOs
                        L7: Cs-(Cds-Cds)(Cds-Cds)CtOs
                        L7: Cs-(Cds-Cdd)(Cds-Cds)CtOs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cds)CtOs
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cds)CtOs
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)CtOs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CtOs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CtOs
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtOs
                L5: Cs-CbCdsCdsOs
                    L6: Cs-(Cds-Od)(Cds-Od)CbOs
                    L6: Cs-(Cds-Od)(Cds-Cd)CbOs
                        L7: Cs-(Cds-Od)(Cds-Cds)CbOs
                        L7: Cs-(Cds-Od)(Cds-Cdd)CbOs
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)CbOs
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)CbOs
                    L6: Cs-(Cds-Cd)(Cds-Cd)CbOs
                        L7: Cs-(Cds-Cds)(Cds-Cds)CbOs
                        L7: Cs-(Cds-Cdd)(Cds-Cds)CbOs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cds)CbOs
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cds)CbOs
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)CbOs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CbOs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CbOs
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbOs
                L5: Cs-CtCtCdsOs
                    L6: Cs-(Cds-Od)CtCtOs
                    L6: Cs-(Cds-Cd)CtCtOs
                        L7: Cs-(Cds-Cds)CtCtOs
                        L7: Cs-(Cds-Cdd)CtCtOs
                            L8: Cs-(Cds-Cdd-Od)CtCtOs
                            L8: Cs-(Cds-Cdd-Cd)CtCtOs
                L5: Cs-CbCtCdsOs
                    L6: Cs-(Cds-Od)CbCtOs
                    L6: Cs-(Cds-Cd)CbCtOs
                        L7: Cs-(Cds-Cds)CbCtOs
                        L7: Cs-(Cds-Cdd)CbCtOs
                            L8: Cs-(Cds-Cdd-Od)CbCtOs
                            L8: Cs-(Cds-Cdd-Cd)CbCtOs
                L5: Cs-CbCbCdsOs
                    L6: Cs-(Cds-Od)CbCbOs
                    L6: Cs-(Cds-Cd)CbCbOs
                        L7: Cs-(Cds-Cds)CbCbOs
                        L7: Cs-(Cds-Cdd)CbCbOs
                            L8: Cs-(Cds-Cdd-Od)CbCbOs
                            L8: Cs-(Cds-Cdd-Cd)CbCbOs
                L5: Cs-CtCtCtOs
                L5: Cs-CbCtCtOs
                L5: Cs-CbCbCtOs
                L5: Cs-CbCbCbOs
            L4: Cs-CCOsOs
                L5: Cs-CsCsOsOs
                L5: Cs-CdsCsOsOs
                    L6: Cs-(Cds-Od)CsOsOs
                    L6: Cs-(Cds-Cd)CsOsOs
                        L7: Cs-(Cds-Cds)CsOsOs
                        L7: Cs-(Cds-Cdd)CsOsOs
                            L8: Cs-(Cds-Cdd-Od)CsOsOs
                            L8: Cs-(Cds-Cdd-Cd)CsOsOs
                L5: Cs-CdsCdsOsOs
                    L6: Cs-(Cds-Od)(Cds-Od)OsOs
                    L6: Cs-(Cds-Od)(Cds-Cd)OsOs
                        L7: Cs-(Cds-Od)(Cds-Cds)OsOs
                        L7: Cs-(Cds-Od)(Cds-Cdd)OsOs
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)OsOs
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)OsOs
                    L6: Cs-(Cds-Cd)(Cds-Cd)OsOs
                        L7: Cs-(Cds-Cds)(Cds-Cds)OsOs
                        L7: Cs-(Cds-Cdd)(Cds-Cds)OsOs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cds)OsOs
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cds)OsOs
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)OsOs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)OsOs
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)OsOs
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)OsOs
                L5: Cs-CtCsOsOs
                L5: Cs-CtCdsOsOs
                    L6: Cs-(Cds-Od)CtOsOs
                    L6: Cs-(Cds-Cd)CtOsOs
                        L7: Cs-(Cds-Cds)CtOsOs
                        L7: Cs-(Cds-Cdd)CtOsOs
                            L8: Cs-(Cds-Cdd-Od)CtOsOs
                            L8: Cs-(Cds-Cdd-Cd)CtOsOs
                L5: Cs-CtCtOsOs
                L5: Cs-CbCsOsOs
                L5: Cs-CbCdsOsOs
                    L6: Cs-(Cds-Od)CbOsOs
                    L6: Cs-(Cds-Cd)CbOsOs
                        L7: Cs-(Cds-Cds)CbOsOs
                        L7: Cs-(Cds-Cdd)CbOsOs
                            L8: Cs-(Cds-Cdd-Od)CbOsOs
                            L8: Cs-(Cds-Cdd-Cd)CbOsOs
                L5: Cs-CbCtOsOs
                L5: Cs-CbCbOsOs
            L4: Cs-COsOsOs
                L5: Cs-CsOsOsOs
                L5: Cs-CdsOsOsOs
                    L6: Cs-(Cds-Od)OsOsOs
                    L6: Cs-(Cds-Cd)OsOsOs
                        L7: Cs-(Cds-Cds)OsOsOs
                        L7: Cs-(Cds-Cdd)OsOsOs
                            L8: Cs-(Cds-Cdd-Od)OsOsOs
                            L8: Cs-(Cds-Cdd-Cd)OsOsOs
                L5: Cs-CtOsOsOs
                L5: Cs-CbOsOsOs
            L4: Cs-OsOsOsOs
            L4: Cs-COsOsH
                L5: Cs-CsOsOsH
                L5: Cs-CdsOsOsH
                    L6: Cs-(Cds-Od)OsOsH
                    L6: Cs-(Cds-Cd)OsOsH
                        L7: Cs-(Cds-Cds)OsOsH
                        L7: Cs-(Cds-Cdd)OsOsH
                            L8: Cs-(Cds-Cdd-Od)OsOsH
                            L8: Cs-(Cds-Cdd-Cd)OsOsH
                L5: Cs-CtOsOsH
                L5: Cs-CbOsOsH
            L4: Cs-CCOsH
                L5: Cs-CsCsOsH
                L5: Cs-CdsCsOsH
                    L6: Cs-(Cds-Od)CsOsH
                    L6: Cs-(Cds-Cd)CsOsH
                        L7: Cs-(Cds-Cds)CsOsH
                        L7: Cs-(Cds-Cdd)CsOsH
                            L8: Cs-(Cds-Cdd-Od)CsOsH
                            L8: Cs-(Cds-Cdd-Cd)CsOsH
                L5: Cs-CdsCdsOsH
                    L6: Cs-(Cds-Od)(Cds-Od)OsH
                    L6: Cs-(Cds-Od)(Cds-Cd)OsH
                        L7: Cs-(Cds-Od)(Cds-Cds)OsH
                        L7: Cs-(Cds-Od)(Cds-Cdd)OsH
                            L8: Cs-(Cds-Od)(Cds-Cdd-Od)OsH
                            L8: Cs-(Cds-Od)(Cds-Cdd-Cd)OsH
                    L6: Cs-(Cds-Cd)(Cds-Cd)OsH
                        L7: Cs-(Cds-Cds)(Cds-Cds)OsH
                        L7: Cs-(Cds-Cdd)(Cds-Cds)OsH
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cds)OsH
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cds)OsH
                        L7: Cs-(Cds-Cdd)(Cds-Cdd)OsH
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)OsH
                            L8: Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)OsH
                            L8: Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)OsH
                L5: Cs-CtCsOsH
                L5: Cs-CtCdsOsH
                    L6: Cs-(Cds-Od)CtOsH
                    L6: Cs-(Cds-Cd)CtOsH
                        L7: Cs-(Cds-Cds)CtOsH
                        L7: Cs-(Cds-Cdd)CtOsH
                            L8: Cs-(Cds-Cdd-Od)CtOsH
                            L8: Cs-(Cds-Cdd-Cd)CtOsH
                L5: Cs-CtCtOsH
                L5: Cs-CbCsOsH
                L5: Cs-CbCdsOsH
                    L6: Cs-(Cds-Od)CbOsH
                    L6: Cs-(Cds-Cd)CbOsH
                        L7: Cs-(Cds-Cds)CbOsH
                        L7: Cs-(Cds-Cdd)CbOsH
                            L8: Cs-(Cds-Cdd-Od)CbOsH
                            L8: Cs-(Cds-Cdd-Cd)CbOsH
                L5: Cs-CbCtOsH
                L5: Cs-CbCbOsH
            L4: Cs-COsHH
                L5: Cs-CsOsHH
                L5: Cs-CdsOsHH
                    L6: Cs-(Cds-Od)OsHH
                    L6: Cs-(Cds-Cd)OsHH
                        L7: Cs-(Cds-Cds)OsHH
                        L7: Cs-(Cds-Cdd)OsHH
                            L8: Cs-(Cds-Cdd-Od)OsHH
                            L8: Cs-(Cds-Cdd-Cd)OsHH
                L5: Cs-CtOsHH
                L5: Cs-CbOsHH
    L2: O
        L3: Od
            L4: Od-Cd
            L4: Od-Od
        L3: Os
            L4: Os-HH
            L4: Os-OsH
            L4: Os-OsOs
            L4: Os-CH
                L5: Os-CtH
                L5: Os-CdsH
                    L6: Os-(Cds-Od)H
                    L6: Os-(Cds-Cd)H
                L5: Os-CsH
                L5: Os-CbH
            L4: Os-OsC
                L5: Os-OsCt
                L5: Os-OsCds
                    L6: Os-Os(Cds-Od)
                    L6: Os-Os(Cds-Cd)
                L5: Os-OsCs
                L5: Os-OsCb
            L4: Os-CC
                L5: Os-CtCt
                L5: Os-CtCds
                    L6: Os-Ct(Cds-Od)
                    L6: Os-Ct(Cds-Cd)
                L5: Os-CtCs
                L5: Os-CtCb
                L5: Os-CdsCds
                    L6: Os-(Cds-Od)(Cds-Od)
                    L6: Os-(Cds-Od)(Cds-Cd)
                    L6: Os-(Cds-Cd)(Cds-Cd)
                L5: Os-CdsCs
                    L6: Os-Cs(Cds-Od)
                    L6: Os-Cs(Cds-Cd)
                L5: Os-CdsCb
                    L6: Os-Cb(Cds-Od)
                    L6: Os-Cb(Cds-Cd)
                L5: Os-CsCs
                L5: Os-CsCb
                L5: Os-CbCb
    L2: Si
    L2: S
"""
)

thermo(
    label="R",
    group=
        """
        1  *  R 0
        """,
    node="C",
    index=0,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="C",
    group=
        """
        1  *  C 0
        """,
    node="Cs-CsCsCsCs",
    index=1,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cbf",
    group=
        """
        1  *  Cbf 0
        """,
    node="Cbf-CbCbCbf",
    index=2,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cbf-CbCbCbf",
    group=
        """
        1  *  Cbf 0 {2,B} {3,B} {4,B}
        2     Cb 0 {1,B}
        3     Cb 0 {1,B}
        4     Cbf 0 {1,B}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.01,3.68,4.2,4.61,5.2,5.7,6.2],"cal/(mol*K)"),
        H298=(4.8,"kcal/mol"),
        S298=(-5,"cal/(mol*K)"),
    ),
    index=3,
    short_comment="Cbf-CbfCbCb STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cbf-CbCbfCbf",
    group=
        """
        1  *  Cbf 0 {2,B} {3,B} {4,B}
        2     Cb 0 {1,B}
        3     Cbf 0 {1,B}
        4     Cbf 0 {1,B}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.01,3.68,4.2,4.61,5.2,5.7,6.2],"cal/(mol*K)"),
        H298=(3.7,"kcal/mol"),
        S298=(-5,"cal/(mol*K)"),
    ),
    index=4,
    short_comment="Cbf-CbfCbfCb STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cbf-CbfCbfCbf",
    group=
        """
        1  *  Cbf 0 {2,B} {3,B} {4,B}
        2     Cbf 0 {1,B}
        3     Cbf 0 {1,B}
        4     Cbf 0 {1,B}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([2,3.11,3.9,4.42,5,5.3,5.7],"cal/(mol*K)"),
        H298=(1.5,"kcal/mol"),
        S298=(1.8,"cal/(mol*K)"),
    ),
    index=5,
    short_comment="Cbf-CbfCbfCbf STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cb",
    group=
        """
        1  *  Cb 0
        """,
    node="Cb-Cs",
    index=6,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cb-H",
    group=
        """
        1  *  Cb 0 {2,S}
        2     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.24,4.44,5.46,6.3,7.54,8.41,9.73],"cal/(mol*K)"),
        H298=(3.3,"kcal/mol"),
        S298=(11.53,"cal/(mol*K)"),
    ),
    index=7,
    short_comment="Cb-H BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cb-Os",
    group=
        """
        1  *  Cb 0 {2,S}
        2     Os 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.9,5.3,6.2,6.6,6.9,6.9,7.07],"cal/(mol*K)"),
        H298=(-0.9,"kcal/mol"),
        S298=(-10.2,"cal/(mol*K)"),
    ),
    index=8,
    short_comment="Cb-O BENSON Cp1500=3D Cp1000*(Cp1500/Cp1000: Cb/Cd)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cb-C",
    group=
        """
        1  *  Cb 0 {2,S}
        2     C 0 {1,S}
        """,
    node="Cb-Cs",
    index=9,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cb-Cs",
    group=
        """
        1  *  Cb 0 {2,S}
        2     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([2.67,3.14,3.68,4.15,4.96,5.44,5.98],"cal/(mol*K)"),
        H298=(5.51,"kcal/mol"),
        S298=(-7.69,"cal/(mol*K)"),
    ),
    index=10,
    short_comment="Cb-Cs BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cb-Cds",
    group=
        """
        1  *  Cb 0 {2,S}
        2     {Cd,CO} 0 {1,S}
        """,
    node="Cb-(Cds-Cds)",
    index=11,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cb-(Cds-Od)",
    group=
        """
        1  *  Cb 0 {2,S}
        2     CO 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.59,3.97,4.38,4.72,5.28,5.61,5.75],"cal/(mol*K)"),
        H298=(3.69,"kcal/mol"),
        S298=(-7.8,"cal/(mol*K)"),
    ),
    index=12,
    short_comment="Enthalpy from Cb-CO, entropies and heat capacities assigned value of Cb-Cd",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cb-(Cds-Cd)",
    group=
        """
        1  *  Cb 0 {2,S}
        2     Cd 0 {1,S}
        """,
    node="Cb-(Cds-Cds)",
    index=13,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cb-(Cds-Cds)",
    group=
        """
        1  *  Cb 0 {2,S}
        2     Cd 0 {1,S} {3,D}
        3     Cd 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.59,3.97,4.38,4.72,5.28,5.61,5.75],"cal/(mol*K)"),
        H298=(5.69,"kcal/mol"),
        S298=(-7.8,"cal/(mol*K)"),
    ),
    index=14,
    short_comment="Cb-Cd STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cb-(Cds-Cdd)",
    group=
        """
        1  *  Cb 0 {2,S}
        2     Cd 0 {1,S} {3,D}
        3     Cdd 0 {2,D}
        """,
    node="Cb-(Cds-Cdd-Cd)",
    index=15,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cb-(Cds-Cdd-Od)",
    group=
        """
        1  *  Cb 0 {2,S}
        2     Cd 0 {1,S} {3,D}
        3     Cdd 0 {2,D} {4,D}
        4     Od 0 {3,D}
        """,
    node="Cb-(Cds-Cds)",
    index=16,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cb-(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cb 0 {2,S}
        2     Cd 0 {1,S} {3,D}
        3     Cdd 0 {2,D} {4,D}
        4     C 0 {3,D}
        """,
    node="Cb-(Cds-Cds)",
    index=17,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cb-Ct",
    group=
        """
        1  *  Cb 0 {2,S}
        2     Ct 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.59,3.97,4.38,4.72,5.28,5.61,5.75],"cal/(mol*K)"),
        H298=(5.69,"kcal/mol"),
        S298=(-7.8,"cal/(mol*K)"),
    ),
    index=18,
    short_comment="Cb-Ct STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cb-Cb",
    group=
        """
        1  *  Cb 0 {2,S}
        2     Cb 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.33,4.22,4.89,5.27,5.76,5.95,6.05],"cal/(mol*K)"),
        H298=(4.96,"kcal/mol"),
        S298=(-8.64,"cal/(mol*K)"),
    ),
    index=19,
    short_comment="Cb-Cb BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Ct",
    group=
        """
        1  *  Ct 0
        """,
    node="Ct-Cs",
    index=20,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Ct-H",
    group=
        """
        1  *  Ct 0 {2,S}
        2     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.28,5.99,6.49,6.87,7.47,7.96,8.85],"cal/(mol*K)"),
        H298=(26.93,"kcal/mol"),
        S298=(24.7,"cal/(mol*K)"),
    ),
    index=21,
    short_comment="Ct-H STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Ct-Os",
    group=
        """
        1  *  Ct 0 {2,S}
        2     Os 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.64,4.39,4.85,5.63,5.66,5.73,5.73],"cal/(mol*K)"),
        H298=(31.4,"kcal/mol"),
        S298=(4.91,"cal/(mol*K)"),
    ),
    index=22,
    short_comment="Ct-O MELIUS / hc#coh !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Ct-C",
    group=
        """
        1  *  Ct 0 {2,S}
        2     C 0 {1,S}
        """,
    node="Ct-Cs",
    index=23,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Ct-Cs",
    group=
        """
        1  *  Ct 0 {2,S}
        2     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.13,3.48,3.81,4.09,4.6,4.92,6.35],"cal/(mol*K)"),
        H298=(27.55,"kcal/mol"),
        S298=(6.35,"cal/(mol*K)"),
    ),
    index=24,
    short_comment="Ct-Cs STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Ct-Cds",
    group=
        """
        1  *  Ct 0 {2,S}
        2     {Cd,CO} 0 {1,S}
        """,
    node="Ct-(Cds-Cds)",
    index=25,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Ct-(Cds-Od)",
    group=
        """
        1  *  Ct 0 {2,S}
        2     CO 0 {1,S}
        """,
    node="Ct-Cs",
    index=26,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Ct-(Cds-Cd)",
    group=
        """
        1  *  Ct 0 {2,S}
        2     Cd 0 {1,S}
        """,
    node="Ct-(Cds-Cds)",
    index=27,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Ct-(Cds-Cds)",
    group=
        """
        1  *  Ct 0 {2,S}
        2     Cd 0 {1,S} {3,D}
        3     Cd 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([2.57,3.54,3.5,4.92,5.34,5.5,5.8],"cal/(mol*K)"),
        H298=(28.2,"kcal/mol"),
        S298=(6.43,"cal/(mol*K)"),
    ),
    index=28,
    short_comment="Ct-Cd STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Ct-(Cds-Cdd)",
    group=
        """
        1  *  Ct 0 {2,S}
        2     Cd 0 {1,S} {3,D}
        3     Cdd 0 {2,D}
        """,
    node="Ct-(Cds-Cdd-Cd)",
    index=29,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Ct-(Cds-Cdd-Od)",
    group=
        """
        1  *  Ct 0 {2,S}
        2     Cd 0 {1,S} {3,D}
        3     Cdd 0 {2,D} {4,D}
        4     Od 0 {3,D}
        """,
    node="Ct-(Cds-Cds)",
    index=30,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Ct-(Cds-Cdd-Cd)",
    group=
        """
        1  *  Ct 0 {2,S}
        2     Cd 0 {1,S} {3,D}
        3     Cdd 0 {2,D} {4,D}
        4     C 0 {3,D}
        """,
    node="Ct-(Cds-Cds)",
    index=31,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Ct-Ct",
    group=
        """
        1  *  Ct 0 {2,S}
        2     Ct 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.54,4.06,4.4,4.64,5,5.23,5.57],"cal/(mol*K)"),
        H298=(25.6,"kcal/mol"),
        S298=(5.88,"cal/(mol*K)"),
    ),
    index=32,
    short_comment="Ct-Ct STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Ct-Cb",
    group=
        """
        1  *  Ct 0 {2,S}
        2     Cb 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([2.57,3.54,4.5,4.92,5.34,5.5,5.8],"cal/(mol*K)"),
        H298=(24.67,"kcal/mol"),
        S298=(6.43,"cal/(mol*K)"),
    ),
    index=33,
    short_comment="Ct-Cb STEIN and FAHR; J. PHYS. CHEM. 1985, 89, 17, 3714",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cdd",
    group=
        """
        1  *  Cdd 0
        """,
    node="Cdd-CdsCds",
    index=34,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cdd-OdOd",
    group=
        """
        1  *  Cdd 0 {2,D} {3,D}
        2     Od 0 {1,D}
        3     Od 0 {1,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([8.91,9.86,10.65,11.31,12.32,12.99,13.93],"cal/(mol*K)"),
        H298=(-94.05,"kcal/mol"),
        S298=(52.46,"cal/(mol*K)"),
    ),
    index=35,
    short_comment="CHEMKIN DATABASE: S(group) = S(CO2) + Rln(2)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cdd-CdOd",
    group=
        """
        1  *  Cdd 0 {2,D} {3,D}
        2     C 0 {1,D}
        3     Od 0 {1,D}
        """,
    node="Cdd-CdsOd",
    index=36,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cdd-CdsOd",
    group=
        """
        1  *  Cdd 0 {2,D} {3,D}
        2     Cd 0 {1,D}
        3     Od 0 {1,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0,0,0,0,0,0,0],"cal/(mol*K)"),
        H298=(0,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=37,
    short_comment="O=C*=C< currently treat the adjacent C as Ck",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cdd-CddOd",
    group=
        """
        1  *  Cdd 0 {2,D} {3,D}
        2     Cdd 0 {1,D}
        3     Od 0 {1,D}
        """,
    node="Cdd-(Cdd-Cd)Od",
    index=38,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cdd-(Cdd-Od)Od",
    group=
        """
        1  *  Cdd 0 {2,D} {3,D}
        2     Cdd 0 {1,D} {4,D}
        3     Od 0 {1,D}
        4     Od 0 {2,D}
        """,
    node="Cdd-CdsOd",
    index=40,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cdd-(Cdd-Cd)Od",
    group=
        """
        1  *  Cdd 0 {2,D} {3,D}
        2     Cdd 0 {1,D} {4,D}
        3     Od 0 {1,D}
        4     C 0 {2,D}
        """,
    node="Cdd-CdsOd",
    index=39,
    short_comment="O=C*=C= currently not defined. Assigned same value as Ca",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cdd-CdCd",
    group=
        """
        1  *  Cdd 0 {2,D} {3,D}
        2     C 0 {1,D}
        3     C 0 {1,D}
        """,
    node="Cdd-CdsCds",
    index=41,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cdd-CddCdd",
    group=
        """
        1  *  Cdd 0 {2,D} {3,D}
        2     Cdd 0 {1,D}
        3     Cdd 0 {1,D}
        """,
    node="Cdd-(Cdd-Cd)(Cdd-Cd)",
    index=42,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cdd-(Cdd-Od)(Cdd-Od)",
    group=
        """
        1  *  Cdd 0 {2,D} {3,D}
        2     Cdd 0 {1,D} {4,D}
        3     Cdd 0 {1,D} {5,D}
        4     Od 0 {2,D}
        5     Od 0 {3,D}
        """,
    node="Cdd-CdsCds",
    index=43,
    short_comment="O=C=C*=C=O, currently not defined. Assigned same value as Ca",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cdd-(Cdd-Od)(Cdd-Cd)",
    group=
        """
        1  *  Cdd 0 {2,D} {3,D}
        2     Cdd 0 {1,D} {4,D}
        3     Cdd 0 {1,D} {5,D}
        4     Od 0 {2,D}
        5     C 0 {3,D}
        """,
    node="Cdd-(Cdd-Od)Cds",
    index=44,
    short_comment="O=C=C*=C=C, currently not defined. Assigned same value as Ca",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cdd-(Cdd-Cd)(Cdd-Cd)",
    group=
        """
        1  *  Cdd 0 {2,D} {3,D}
        2     Cdd 0 {1,D} {4,D}
        3     Cdd 0 {1,D} {5,D}
        4     C 0 {2,D}
        5     C 0 {3,D}
        """,
    node="Cdd-CdsCds",
    index=45,
    short_comment="C=C=C*=C=C, currently not defined. Assigned same value as Ca",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cdd-CddCds",
    group=
        """
        1  *  Cdd 0 {2,D} {3,D}
        2     Cdd 0 {1,D}
        3     Cd 0 {1,D}
        """,
    node="Cdd-(Cdd-Cd)(Cdd-Cd)",
    index=46,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cdd-(Cdd-Od)Cds",
    group=
        """
        1  *  Cdd 0 {2,D} {3,D}
        2     Cdd 0 {1,D} {4,D}
        3     Cd 0 {1,D}
        4     Od 0 {2,D}
        """,
    node="Cdd-CdsCds",
    index=47,
    short_comment="O=C=C*=C<, currently not defined. Assigned same value as Ca",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cdd-(Cdd-Cd)Cds",
    group=
        """
        1  *  Cdd 0 {2,D} {3,D}
        2     Cdd 0 {1,D} {4,D}
        3     Cd 0 {1,D}
        4     C 0 {2,D}
        """,
    node="Cdd-CdsCds",
    index=48,
    short_comment="C=C=C*=C<, currently not defined. Assigned same value as Ca",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cdd-CdsCds",
    group=
        """
        1  *  Cdd 0 {2,D} {3,D}
        2     Cd 0 {1,D}
        3     Cd 0 {1,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.9,4.4,4.7,5,5.3,5.5,5.7],"cal/(mol*K)"),
        H298=(34.2,"kcal/mol"),
        S298=(6,"cal/(mol*K)"),
    ),
    index=49,
    short_comment="Benson's Ca",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds",
    group=
        """
        1  *  {Cd,CO} 0
        """,
    node="Cds-CdsCsCs",
    index=50,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdHH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     H 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([8.47,9.38,10.46,11.52,13.37,14.81,14.81],"cal/(mol*K)"),
        H298=(-25.95,"kcal/mol"),
        S298=(53.68,"cal/(mol*K)"),
    ),
    index=51,
    short_comment="CO-HH BENSON !!!WARNING! Cp1500 value taken as Cp1000, S(group) = S(CH2O) + Rln(2)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdOsH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Os 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([7.03,7.87,8.82,9.68,11.2,12.2,12.2],"cal/(mol*K)"),
        H298=(-32.1,"kcal/mol"),
        S298=(34.9,"cal/(mol*K)"),
    ),
    index=52,
    short_comment="CO-OH BENSON !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdOsOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.97,6.7,7.4,8.02,8.87,9.36,9.36],"cal/(mol*K)"),
        H298=(-31.45,"kcal/mol"),
        S298=(10.78,"cal/(mol*K)"),
    ),
    index=53,
    short_comment="CO-OO BOZZELLI 8/91, S CO/C/O, Hf PEDLEY ccoc*oocc Bsn Hf-24 !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     C 0 {1,S}
        4     H 0 {1,S}
        """,
    node="Cds-OdCsH",
    index=54,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCsH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([7.03,7.87,8.82,9.68,11.2,12.2,12.2],"cal/(mol*K)"),
        H298=(-29.1,"kcal/mol"),
        S298=(34.9,"cal/(mol*K)"),
    ),
    index=55,
    short_comment="CO-CsH BENSON !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCdsH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     {Cd,CO} 0 {1,S}
        4     H 0 {1,S}
        """,
    node="Cds-Od(Cds-Cds)H",
    index=56,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Od)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     CO 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([7.45,8.77,9.82,10.7,12.14,12.9,12.9],"cal/(mol*K)"),
        H298=(-25.3,"kcal/mol"),
        S298=(33.4,"cal/(mol*K)"),
    ),
    index=57,
    short_comment="CO-COH Hf BENSON S,Cp =3D CO/Cs/H-del(Cd/Cd/H-Cd/Cs/H) !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cd)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S}
        4     H 0 {1,S}
        """,
    node="Cds-Od(Cds-Cds)H",
    index=58,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cds)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     H 0 {1,S}
        5     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([7.45,8.77,9.82,10.7,12.14,12.9,12.9],"cal/(mol*K)"),
        H298=(-30.9,"kcal/mol"),
        S298=(33.4,"cal/(mol*K)"),
    ),
    index=59,
    short_comment="CO-CdH Hf BOZZELLI S,Cp =3D CO/C/H-del(cd syst) !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cdd)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     H 0 {1,S}
        5     Cdd 0 {3,D}
        """,
    node="Cds-Od(Cds-Cdd-Cd)H",
    index=60,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cdd-Od)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     H 0 {1,S}
        5     Cdd 0 {3,D} {6,D}
        6     Od 0 {5,D}
        """,
    node="Cds-Od(Cds-Cds)H",
    index=61,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cdd-Cd)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     H 0 {1,S}
        5     Cdd 0 {3,D} {6,D}
        6     C 0 {5,D}
        """,
    node="Cds-Od(Cds-Cds)H",
    index=62,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCtH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Ct 0 {1,S}
        4     H 0 {1,S}
        """,
    node="Cds-Od(Cds-Cds)H",
    index=63,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCbH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cb 0 {1,S}
        4     H 0 {1,S}
        """,
    node="Cds-Od(Cds-Cds)H",
    index=64,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     C 0 {1,S}
        4     Os 0 {1,S}
        """,
    node="Cds-OdCsOs",
    index=65,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCsOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([6.1,6.7,7.4,8.02,8.87,9.36,9.36],"cal/(mol*K)"),
        H298=(-35.1,"kcal/mol"),
        S298=(10.04,"cal/(mol*K)"),
    ),
    index=66,
    short_comment="CO-OCs Hf BENSON S STULL !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCdsOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     {Cd,CO} 0 {1,S}
        4     Os 0 {1,S}
        """,
    node="Cds-Od(Cds-Cds)Os",
    index=67,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Od)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     CO 0 {1,S}
        4     Os 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.46,6.32,7.17,7.88,9,9.77,9.77],"cal/(mol*K)"),
        H298=(-29.3,"kcal/mol"),
        S298=(14.6,"cal/(mol*K)"),
    ),
    index=68,
    short_comment="CO-OCO Hf,S BOZZELLI Cp BENSON !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cd)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S}
        4     Os 0 {1,S}
        """,
    node="Cds-Od(Cds-Cds)Os",
    index=69,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cds)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Os 0 {1,S}
        5     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.97,6.7,7.4,8.02,8.87,9.36,9.36],"cal/(mol*K)"),
        H298=(-32.1,"kcal/mol"),
        S298=(14.78,"cal/(mol*K)"),
    ),
    index=70,
    short_comment="CO-OCd RPS + S Coreected !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cdd)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Os 0 {1,S}
        5     Cdd 0 {3,D}
        """,
    node="Cds-Od(Cds-Cdd-Cd)Os",
    index=71,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cdd-Od)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Os 0 {1,S}
        5     Cdd 0 {3,D} {6,D}
        6     Od 0 {5,D}
        """,
    node="Cds-Od(Cds-Cds)Os",
    index=72,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cdd-Cd)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Os 0 {1,S}
        5     Cdd 0 {3,D} {6,D}
        6     C 0 {5,D}
        """,
    node="Cds-Od(Cds-Cds)Os",
    index=73,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCtOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        """,
    node="Cds-Od(Cds-Cds)Os",
    index=74,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCbOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cb 0 {1,S}
        4     Os 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.97,6.7,7.4,8.02,8.87,9.36,9.36],"cal/(mol*K)"),
        H298=(-36.6,"kcal/mol"),
        S298=(14.78,"cal/(mol*K)"),
    ),
    index=75,
    short_comment="CO-OCb RPS + S Coreected !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCC",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     C 0 {1,S}
        4     C 0 {1,S}
        """,
    node="Cds-OdCsCs",
    index=76,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCsCs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.59,6.32,7.09,7.76,8.89,9.61,9.61],"cal/(mol*K)"),
        H298=(-31.4,"kcal/mol"),
        S298=(15.01,"cal/(mol*K)"),
    ),
    index=77,
    short_comment="CO-CsCs BENSON !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCdsCs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     {Cd,CO} 0 {1,S}
        4     Cs 0 {1,S}
        """,
    node="Cds-Od(Cds-Cds)Cs",
    index=78,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Od)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     CO 0 {1,S}
        4     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.46,6.32,7.17,7.88,9,9.77,9.77],"cal/(mol*K)"),
        H298=(-29.1,"kcal/mol"),
        S298=(14.6,"cal/(mol*K)"),
    ),
    index=79,
    short_comment="CO-COCs Hf,S BOZZELLI Cp BENSON !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cd)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S}
        4     Cs 0 {1,S}
        """,
    node="Cds-Od(Cds-Cds)Cs",
    index=80,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cds)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cs 0 {1,S}
        5     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.46,6.32,7.17,7.88,9,9.77,9.77],"cal/(mol*K)"),
        H298=(-30.9,"kcal/mol"),
        S298=(14.6,"cal/(mol*K)"),
    ),
    index=81,
    short_comment="CO-CdCs Hf BENSON =3D CO/Cb/C S,Cp !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cdd)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cs 0 {1,S}
        5     Cdd 0 {3,D}
        """,
    node="Cds-Od(Cds-Cdd-Cd)Cs",
    index=82,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cdd-Od)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cs 0 {1,S}
        5     Cdd 0 {3,D} {6,D}
        6     Od 0 {5,D}
        """,
    node="Cds-Od(Cds-Cds)Cs",
    index=83,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cdd-Cd)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cs 0 {1,S}
        5     Cdd 0 {3,D} {6,D}
        6     C 0 {5,D}
        """,
    node="Cds-Od(Cds-Cds)Cs",
    index=84,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCdsCds",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     {Cd,CO} 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        """,
    node="Cds-Od(Cds-Cds)(Cds-Cds)",
    index=85,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Od)(Cds-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     CO 0 {1,S}
        4     CO 0 {1,S}
        """,
    node="Cds-OdCsCs",
    index=86,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cd)(Cds-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S}
        4     CO 0 {1,S}
        """,
    node="Cds-Od(Cds-Cds)(Cds-Od)",
    index=87,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cds)(Cds-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     CO 0 {1,S}
        5     Cd 0 {3,D}
        """,
    node="Cds-Od(Cds-Od)Cs",
    index=88,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cdd)(Cds-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     CO 0 {1,S}
        5     Cdd 0 {3,D}
        """,
    node="Cds-Od(Cds-Cdd-Cd)(Cds-Od)",
    index=89,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cdd-Od)(Cds-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     CO 0 {1,S}
        5     Cdd 0 {3,D} {6,D}
        6     Od 0 {5,D}
        """,
    node="Cds-Od(Cds-Cdd-Od)Cs",
    index=90,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cdd-Cd)(Cds-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     CO 0 {1,S}
        5     Cdd 0 {3,D} {6,D}
        6     C 0 {5,D}
        """,
    node="Cds-Od(Cds-Cds)(Cds-Od)",
    index=91,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cd)(Cds-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S}
        4     Cd 0 {1,S}
        """,
    node="Cds-Od(Cds-Cds)(Cds-Cds)",
    index=92,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cds)(Cds-Cds)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cd 0 {1,S} {6,D}
        5     Cd 0 {3,D}
        6     Cd 0 {4,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.46,6.32,7.17,7.88,9,9.77,9.77],"cal/(mol*K)"),
        H298=(-30.9,"kcal/mol"),
        S298=(14.6,"cal/(mol*K)"),
    ),
    index=93,
    short_comment="CO-CdCd Estimate BOZZELLI !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cdd)(Cds-Cds)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cd 0 {1,S} {6,D}
        5     Cdd 0 {3,D}
        6     Cd 0 {4,D}
        """,
    node="Cds-Od(Cds-Cdd-Cd)(Cds-Cds)",
    index=94,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cdd-Od)(Cds-Cds)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cd 0 {1,S} {6,D}
        5     Cdd 0 {3,D} {7,D}
        6     Cd 0 {4,D}
        7     Od 0 {5,D}
        """,
    node="Cds-Od(Cds-Cdd-Od)Cs",
    index=95,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cdd-Cd)(Cds-Cds)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cd 0 {1,S} {6,D}
        5     Cdd 0 {3,D} {7,D}
        6     Cd 0 {4,D}
        7     C 0 {5,D}
        """,
    node="Cds-Od(Cds-Cds)(Cds-Cds)",
    index=96,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cdd)(Cds-Cdd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cd 0 {1,S} {6,D}
        5     Cdd 0 {3,D}
        6     Cdd 0 {4,D}
        """,
    node="Cds-Od(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    index=97,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cdd-Od)(Cds-Cdd-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cd 0 {1,S} {6,D}
        5     Cdd 0 {3,D} {7,D}
        6     Cdd 0 {4,D} {8,D}
        7     Od 0 {5,D}
        8     Od 0 {6,D}
        """,
    node="Cds-Od(Cds-Cds)(Cds-Cds)",
    index=98,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cdd-Cd)(Cds-Cdd-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cd 0 {1,S} {6,D}
        5     Cdd 0 {3,D} {7,D}
        6     Cdd 0 {4,D} {8,D}
        7     C 0 {5,D}
        8     Od 0 {6,D}
        """,
    node="Cds-Od(Cds-Cdd-Od)(Cds-Cds)",
    index=99,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Od(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cd 0 {1,S} {6,D}
        5     Cdd 0 {3,D} {7,D}
        6     Cdd 0 {4,D} {8,D}
        7     C 0 {5,D}
        8     C 0 {6,D}
        """,
    node="Cds-Od(Cds-Cds)(Cds-Cds)",
    index=100,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCtCs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        """,
    node="Cds-Od(Cds-Cds)Cs",
    index=101,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCtCds",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Ct 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        """,
    node="Cds-OdCt(Cds-Cds)",
    index=102,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCt(Cds-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Ct 0 {1,S}
        4     CO 0 {1,S}
        """,
    node="Cds-Od(Cds-Cds)(Cds-Od)",
    index=103,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCt(Cds-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Ct 0 {1,S}
        4     Cd 0 {1,S}
        """,
    node="Cds-OdCt(Cds-Cds)",
    index=104,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCt(Cds-Cds)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Ct 0 {1,S}
        4     Cd 0 {1,S} {5,D}
        5     Cd 0 {4,D}
        """,
    node="Cds-Od(Cds-Cds)(Cds-Cds)",
    index=105,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCt(Cds-Cdd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Ct 0 {1,S}
        4     Cd 0 {1,S} {5,D}
        5     Cdd 0 {4,D}
        """,
    node="Cds-OdCt(Cds-Cdd-Cd)",
    index=106,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCt(Cds-Cdd-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Ct 0 {1,S}
        4     Cd 0 {1,S} {5,D}
        5     Cdd 0 {4,D} {6,D}
        6     Od 0 {5,D}
        """,
    node="Cds-Od(Cds-Cdd-Od)(Cds-Cds)",
    index=107,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCt(Cds-Cdd-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Ct 0 {1,S}
        4     Cd 0 {1,S} {5,D}
        5     Cdd 0 {4,D} {6,D}
        6     C 0 {5,D}
        """,
    node="Cds-OdCt(Cds-Cds)",
    index=108,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCtCt",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        """,
    node="Cds-Od(Cds-Cds)(Cds-Cds)",
    index=109,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCbCs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        """,
    node="Cds-Od(Cds-Cds)Cs",
    index=110,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCbCds",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cb 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        """,
    node="Cds-OdCb(Cds-Cds)",
    index=111,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCb(Cds-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cb 0 {1,S}
        4     CO 0 {1,S}
        """,
    node="Cds-Od(Cds-Cds)(Cds-Od)",
    index=112,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCb(Cds-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cb 0 {1,S}
        4     Cd 0 {1,S}
        """,
    node="Cds-OdCb(Cds-Cds)",
    index=113,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCb(Cds-Cds)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cb 0 {1,S}
        4     Cd 0 {1,S} {5,D}
        5     Cd 0 {4,D}
        """,
    node="Cds-Od(Cds-Cds)(Cds-Cds)",
    index=114,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCb(Cds-Cdd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cb 0 {1,S}
        4     Cd 0 {1,S} {5,D}
        5     Cdd 0 {4,D}
        """,
    node="Cds-OdCb(Cds-Cdd-Cd)",
    index=115,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCb(Cds-Cdd-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cb 0 {1,S}
        4     Cd 0 {1,S} {5,D}
        5     Cdd 0 {4,D} {6,D}
        6     Od 0 {5,D}
        """,
    node="Cds-Od(Cds-Cdd-Od)(Cds-Cds)",
    index=116,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCb(Cds-Cdd-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cb 0 {1,S}
        4     Cd 0 {1,S} {5,D}
        5     Cdd 0 {4,D} {6,D}
        6     C 0 {5,D}
        """,
    node="Cds-OdCb(Cds-Cds)",
    index=117,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCbCt",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        """,
    node="Cds-OdCt(Cds-Cds)",
    index=118,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-OdCbCb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Od 0 {1,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        """,
    node="Cds-Od(Cds-Cds)(Cds-Cds)",
    index=119,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdHH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     C 0 {1,D}
        3     H 0 {1,S}
        4     H 0 {1,S}
        """,
    node="Cds-CdsHH",
    index=120,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsHH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     H 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.1,6.36,7.51,8.5,10.07,11.27,13.19],"cal/(mol*K)"),
        H298=(6.26,"kcal/mol"),
        S298=(27.61,"cal/(mol*K)"),
    ),
    index=121,
    short_comment="Cd-HH BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddHH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     H 0 {1,S}
        4     H 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)HH",
    index=122,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)HH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     Od 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([12.08,13.91,15.4,16.64,18.61,20.1,22.47],"cal/(mol*K)"),
        H298=(-11.34,"kcal/mol"),
        S298=(57.47,"cal/(mol*K)"),
    ),
    index=123,
    short_comment="{CCO/H2} RAMAN & GREEN JPCA 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)HH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-CdsHH",
    index=124,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdOsH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     C 0 {1,D}
        3     Os 0 {1,S}
        4     H 0 {1,S}
        """,
    node="Cds-CdsOsH",
    index=125,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsOsH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Os 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.75,6.46,7.64,8.35,9.1,9.56,10.46],"cal/(mol*K)"),
        H298=(2.03,"kcal/mol"),
        S298=(6.2,"cal/(mol*K)"),
    ),
    index=126,
    short_comment="Cd-OH BOZZELLI Hf vin-oh RADOM + C/Cd/H, S&Cp LAY",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddOsH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     Os 0 {1,S}
        4     H 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)OsH",
    index=127,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)OsH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Os 0 {1,S}
        4     H 0 {1,S}
        5     Od 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([11.29,13.67,15.1,16.1,17.36,18.25,19.75],"cal/(mol*K)"),
        H298=(2.11,"kcal/mol"),
        S298=(38.17,"cal/(mol*K)"),
    ),
    index=128,
    short_comment="{CCO/O/H} RAMAN & GREEN JPCA 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)OsH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Os 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-CdsOsH",
    index=129,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdOsOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     C 0 {1,D}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        """,
    node="Cds-CdsOsOs",
    index=130,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsOsOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        """,
    node="Cds-CdsCsCs",
    index=131,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddOsOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)OsOs",
    index=132,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)OsOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     Od 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([11.56,15.58,17.69,18.67,18.78,18.4,18.01],"cal/(mol*K)"),
        H298=(2.403,"kcal/mol"),
        S298=(13.42,"cal/(mol*K)"),
    ),
    index=133,
    short_comment="{CCO/O2} RAMAN & GREEN JPCA 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)OsOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-CdsOsOs",
    index=134,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdCH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     C 0 {1,D}
        3     C 0 {1,S}
        4     H 0 {1,S}
        """,
    node="Cds-CdsCsH",
    index=135,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCsH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.16,5.03,5.81,6.5,7.65,8.45,9.62],"cal/(mol*K)"),
        H298=(8.59,"kcal/mol"),
        S298=(7.97,"cal/(mol*K)"),
    ),
    index=136,
    short_comment="Cd-CsH BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCdsH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     {Cd,CO} 0 {1,S}
        4     H 0 {1,S}
        """,
    node="Cds-Cds(Cds-Cds)H",
    index=137,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Od)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     CO 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.46,5.79,6.75,7.42,8.35,9.11,10.09],"cal/(mol*K)"),
        H298=(4.32,"kcal/mol"),
        S298=(6.38,"cal/(mol*K)"),
    ),
    index=138,
    short_comment="Cd-COH BOZZELLI lit rev Jul91 S,Cp Cd/Cd/H",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cd)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S}
        4     H 0 {1,S}
        """,
    node="Cds-Cds(Cds-Cds)H",
    index=139,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cds)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     H 0 {1,S}
        5     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.46,5.79,6.75,7.42,8.35,8.99,9.98],"cal/(mol*K)"),
        H298=(6.78,"kcal/mol"),
        S298=(6.38,"cal/(mol*K)"),
    ),
    index=140,
    short_comment="Cd-CdH BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cdd)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     H 0 {1,S}
        5     Cdd 0 {3,D}
        """,
    node="Cds-Cds(Cds-Cdd-Cd)H",
    index=141,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cdd-Od)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     H 0 {1,S}
        5     Cdd 0 {3,D} {6,D}
        6     Od 0 {5,D}
        """,
    node="Cds-Cds(Cds-Cds)H",
    index=142,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cdd-Cd)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     H 0 {1,S}
        5     Cdd 0 {3,D} {6,D}
        6     C 0 {5,D}
        """,
    node="Cds-Cds(Cds-Cds)H",
    index=143,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCtH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Ct 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.46,5.79,6.75,7.42,8.35,8.99,9.98],"cal/(mol*K)"),
        H298=(6.78,"kcal/mol"),
        S298=(6.38,"cal/(mol*K)"),
    ),
    index=144,
    short_comment="Cd-CtH BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCbH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cb 0 {1,S}
        4     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.46,5.79,6.75,7.42,8.35,8.99,9.98],"cal/(mol*K)"),
        H298=(6.78,"kcal/mol"),
        S298=(6.38,"cal/(mol*K)"),
    ),
    index=145,
    short_comment="Cd-CbH BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddCsH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)CsH",
    index=146,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)CsH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        5     Od 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([10.31,11.72,12.94,13.98,15.71,16.95,18.78],"cal/(mol*K)"),
        H298=(-4.947,"kcal/mol"),
        S298=(40.04,"cal/(mol*K)"),
    ),
    index=147,
    short_comment="{CCO/H/C} RAMAN & GREEN JPCA 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)CsH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-CdsCsH",
    index=148,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddCdsH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     {Cd,CO} 0 {1,S}
        4     H 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cds)H",
    index=149,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Od)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     CO 0 {1,S}
        4     H 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)CsH",
    index=150,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cd)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S}
        4     H 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)H",
    index=151,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cds)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     H 0 {1,S}
        5     Od 0 {2,D}
        6     Cd 0 {3,D}
        """,
    node="Cds-(Cdd-Od)CsH",
    index=152,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     H 0 {1,S}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cdd-Cd)H",
    index=153,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd-Od)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     H 0 {1,S}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([10.55,12.41,13.82,14.91,16.51,17.62,19.24],"cal/(mol*K)"),
        H298=(-4.998,"kcal/mol"),
        S298=(39.06,"cal/(mol*K)"),
    ),
    index=154,
    short_comment="{CCO/H/CCO} RAMAN & GREEN JPCA 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd-Cd)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     H 0 {1,S}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)H",
    index=155,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Od)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     CO 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-Cds(Cds-Od)H",
    index=156,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cd)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cds)H",
    index=157,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cds)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     H 0 {1,S}
        5     C 0 {2,D}
        6     Cd 0 {3,D}
        """,
    node="Cds-Cds(Cds-Cds)H",
    index=158,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     H 0 {1,S}
        5     C 0 {2,D}
        6     Cdd 0 {3,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cdd-Cd)H",
    index=159,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd-Od)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     H 0 {1,S}
        5     C 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cds-Cds(Cds-Cdd-Od)H",
    index=160,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd-Cd)H",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     H 0 {1,S}
        5     C 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cds)H",
    index=161,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddCtH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     Ct 0 {1,S}
        4     H 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)CtH",
    index=162,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)CtH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Ct 0 {1,S}
        4     H 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)H",
    index=163,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)CtH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Ct 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-CdsCtH",
    index=164,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddCbH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     Cb 0 {1,S}
        4     H 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)CbH",
    index=165,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)CbH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cb 0 {1,S}
        4     H 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)H",
    index=166,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)CbH",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cb 0 {1,S}
        4     H 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-CdsCbH",
    index=167,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdCO",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     C 0 {1,D}
        3     C 0 {1,S}
        4     O 0 {1,S}
        """,
    node="Cds-CdsCsOs",
    index=168,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCsOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.59,4.56,5.04,5.3,5.84,6.07,6.16],"cal/(mol*K)"),
        H298=(3.03,"kcal/mol"),
        S298=(-12.32,"cal/(mol*K)"),
    ),
    index=169,
    short_comment="Cd-OCs BOZZELLI-RADOM vin-oh and del (ccoh-ccohc)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCdsOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     {Cd,CO} 0 {1,S}
        4     Os 0 {1,S}
        """,
    node="Cds-Cds(Cds-Cds)Os",
    index=170,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Od)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     CO 0 {1,S}
        4     Os 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.4,5.37,5.93,6.18,6.5,6.62,6.72],"cal/(mol*K)"),
        H298=(5.13,"kcal/mol"),
        S298=(-14.6,"cal/(mol*K)"),
    ),
    index=171,
    short_comment="Cd-OCO adj BENSON for RADOM c*coh",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cd)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S}
        4     Os 0 {1,S}
        """,
    node="Cds-Cds(Cds-Cds)Os",
    index=172,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cds)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Os 0 {1,S}
        5     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.4,5.37,5.93,6.18,6.5,6.62,6.72],"cal/(mol*K)"),
        H298=(1.5,"kcal/mol"),
        S298=(-14.4,"cal/(mol*K)"),
    ),
    index=173,
    short_comment="Cd-OCd jwb need calc",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cdd)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Os 0 {1,S}
        5     Cdd 0 {3,D}
        """,
    node="Cds-Cds(Cds-Cdd-Cd)Os",
    index=174,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cdd-Od)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Os 0 {1,S}
        5     Cdd 0 {3,D} {6,D}
        6     Od 0 {5,D}
        """,
    node="Cds-Cds(Cds-Cds)Os",
    index=175,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cdd-Cd)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Os 0 {1,S}
        5     Cdd 0 {3,D} {6,D}
        6     C 0 {5,D}
        """,
    node="Cds-Cds(Cds-Cds)Os",
    index=176,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCtOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        """,
    node="Cds-Cds(Cds-Cds)Os",
    index=177,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCbOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cb 0 {1,S}
        4     Os 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.4,5.37,5.93,6.18,6.5,6.62,6.72],"cal/(mol*K)"),
        H298=(1.5,"kcal/mol"),
        S298=(-14.4,"cal/(mol*K)"),
    ),
    index=178,
    short_comment="Cd-OCb jwb need calc",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddCsOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)CsOs",
    index=179,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)CsOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     Od 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([10.91,12.65,13.59,14.22,15,15.48,16.28],"cal/(mol*K)"),
        H298=(3.273,"kcal/mol"),
        S298=(18.58,"cal/(mol*K)"),
    ),
    index=180,
    short_comment="{CCO/O/C} RAMAN & GREEN JPCA 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)CsOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-CdsCsOs",
    index=181,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddCdsOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     {Cd,CO} 0 {1,S}
        4     Os 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cds)Os",
    index=182,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Od)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     CO 0 {1,S}
        4     Os 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)CsOs",
    index=183,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cd)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S}
        4     Os 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)Os",
    index=184,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cds)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Os 0 {1,S}
        5     Od 0 {2,D}
        6     Cd 0 {3,D}
        """,
    node="Cds-(Cdd-Od)CsOs",
    index=185,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Os 0 {1,S}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cdd-Cd)Os",
    index=186,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd-Od)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Os 0 {1,S}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([11.01,12.97,14.17,14.97,15.8,16.26,16.88],"cal/(mol*K)"),
        H298=(1.607,"kcal/mol"),
        S298=(17.73,"cal/(mol*K)"),
    ),
    index=187,
    short_comment="{CCO/O/CCO} RAMAN & GREEN JPCA 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd-Cd)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Os 0 {1,S}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)Os",
    index=188,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cd)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S}
        4     Os 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cds)Os",
    index=189,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cds)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Os 0 {1,S}
        5     C 0 {2,D}
        6     Cd 0 {3,D}
        """,
    node="Cds-Cds(Cds-Cds)Os",
    index=190,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Os 0 {1,S}
        5     C 0 {2,D}
        6     Cdd 0 {3,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cdd-Cd)Os",
    index=191,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd-Od)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Os 0 {1,S}
        5     C 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cds-Cds(Cds-Cdd-Od)Os",
    index=192,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd-Cd)Os",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Os 0 {1,S}
        5     C 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cds)Os",
    index=193,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddCtOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)CtOs",
    index=194,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)CtOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)Os",
    index=195,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)CtOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-CdsCtOs",
    index=196,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddCbOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     Cb 0 {1,S}
        4     Os 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)CbOs",
    index=197,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)CbOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cb 0 {1,S}
        4     Os 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)Os",
    index=198,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)CbOs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cb 0 {1,S}
        4     Os 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-CdsCbOs",
    index=199,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdCC",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     C 0 {1,D}
        3     C 0 {1,S}
        4     C 0 {1,S}
        """,
    node="Cds-CdsCsCs",
    index=200,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCsCs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.1,4.61,4.99,5.26,5.8,6.08,6.36],"cal/(mol*K)"),
        H298=(10.34,"kcal/mol"),
        S298=(-12.7,"cal/(mol*K)"),
    ),
    index=201,
    short_comment="Cd-CsCs BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCdsCs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     {Cd,CO} 0 {1,S}
        4     Cs 0 {1,S}
        """,
    node="Cds-Cds(Cds-Cds)Cs",
    index=202,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Od)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     CO 0 {1,S}
        4     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.4,5.37,5.93,6.18,6.5,6.62,6.72],"cal/(mol*K)"),
        H298=(7.5,"kcal/mol"),
        S298=(-14.6,"cal/(mol*K)"),
    ),
    index=203,
    short_comment="Cd-COCs BENSON Hf, Cd/C/Cd =3D S,Cp",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cd)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S}
        4     Cs 0 {1,S}
        """,
    node="Cds-Cds(Cds-Cds)Cs",
    index=204,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cds)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cs 0 {1,S}
        5     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.4,5.37,5.93,6.18,6.5,6.62,6.72],"cal/(mol*K)"),
        H298=(8.88,"kcal/mol"),
        S298=(-14.6,"cal/(mol*K)"),
    ),
    index=205,
    short_comment="Cd-CdCs BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cdd)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cs 0 {1,S}
        5     Cdd 0 {3,D}
        """,
    node="Cds-Cds(Cds-Cdd-Cd)Cs",
    index=206,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cdd-Od)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cs 0 {1,S}
        5     Cdd 0 {3,D} {6,D}
        6     Od 0 {5,D}
        """,
    node="Cds-Cds(Cds-Cds)Cs",
    index=207,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cdd-Cd)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cs 0 {1,S}
        5     Cdd 0 {3,D} {6,D}
        6     C 0 {5,D}
        """,
    node="Cds-Cds(Cds-Cds)Cs",
    index=208,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCdsCds",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     {Cd,CO} 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        """,
    node="Cds-Cds(Cds-Cds)(Cds-Cds)",
    index=209,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Od)(Cds-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     CO 0 {1,S}
        4     CO 0 {1,S}
        """,
    node="Cds-CdsCsCs",
    index=210,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Od)(Cds-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     CO 0 {1,S}
        4     Cd 0 {1,S}
        """,
    node="Cds-Cds(Cds-Od)(Cds-Cds)",
    index=211,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Od)(Cds-Cds)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {5,D}
        5     Cd 0 {4,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.7,6.13,6.87,7.1,7.2,7.16,7.06],"cal/(mol*K)"),
        H298=(4.6,"kcal/mol"),
        S298=(-16.5,"cal/(mol*K)"),
    ),
    index=212,
    short_comment="Cd-COCd from CD/CD2/ jwb est 6/97",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Od)(Cds-Cdd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {5,D}
        5     Cdd 0 {4,D}
        """,
    node="Cds-Cds(Cds-Od)(Cds-Cdd-Cd)",
    index=213,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Od)(Cds-Cdd-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {5,D}
        5     Cdd 0 {4,D} {6,D}
        6     Od 0 {5,D}
        """,
    node="Cds-Cds(Cds-Cdd-Od)Cs",
    index=214,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Od)(Cds-Cdd-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {5,D}
        5     Cdd 0 {4,D} {6,D}
        6     C 0 {5,D}
        """,
    node="Cds-Cds(Cds-Od)(Cds-Cds)",
    index=215,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cd)(Cds-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S}
        4     Cd 0 {1,S}
        """,
    node="Cds-Cds(Cds-Cds)(Cds-Cds)",
    index=216,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cds)(Cds-Cds)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cd 0 {1,S} {6,D}
        5     Cd 0 {3,D}
        6     Cd 0 {4,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([1.9,2.69,3.5,4.28,5.57,6.21,7.37],"cal/(mol*K)"),
        H298=(4.6,"kcal/mol"),
        S298=(-15.67,"cal/(mol*K)"),
    ),
    index=217,
    short_comment="Cd-CdCd Hf=3D est S,Cp mopac nov99",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cds)(Cds-Cdd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cd 0 {1,S} {6,D}
        5     Cd 0 {3,D}
        6     Cdd 0 {4,D}
        """,
    node="Cds-Cds(Cds-Cds)(Cds-Cdd-Cd)",
    index=218,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cds)(Cds-Cdd-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cd 0 {1,S} {6,D}
        5     Cd 0 {3,D}
        6     Cdd 0 {4,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cds-Cds(Cds-Cdd-Od)Cs",
    index=219,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cds)(Cds-Cdd-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cd 0 {1,S} {6,D}
        5     Cd 0 {3,D}
        6     Cdd 0 {4,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cds-Cds(Cds-Cds)(Cds-Cds)",
    index=220,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cdd)(Cds-Cdd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cd 0 {1,S} {6,D}
        5     Cdd 0 {3,D}
        6     Cdd 0 {4,D}
        """,
    node="Cds-Cds(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    index=221,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cdd-Od)(Cds-Cdd-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cd 0 {1,S} {6,D}
        5     Cdd 0 {3,D} {7,D}
        6     Cdd 0 {4,D} {8,D}
        7     Od 0 {5,D}
        8     Od 0 {6,D}
        """,
    node="Cds-Cds(Cds-Cds)(Cds-Cds)",
    index=222,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cdd-Od)(Cds-Cdd-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cd 0 {1,S} {6,D}
        5     Cdd 0 {3,D} {7,D}
        6     Cdd 0 {4,D} {8,D}
        7     Od 0 {5,D}
        8     C 0 {6,D}
        """,
    node="Cds-Cds(Cds-Cds)(Cds-Cdd-Od)",
    index=223,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cd 0 {1,S} {6,D}
        5     Cdd 0 {3,D} {7,D}
        6     Cdd 0 {4,D} {8,D}
        7     C 0 {5,D}
        8     C 0 {6,D}
        """,
    node="Cds-Cds(Cds-Cds)(Cds-Cdd-Cd)",
    index=224,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCtCs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.5,3.88,4.88,4.18,4.86,5.4,6.01],"cal/(mol*K)"),
        H298=(8.11,"kcal/mol"),
        S298=(-13.02,"cal/(mol*K)"),
    ),
    index=225,
    short_comment="Cd-CtCs RAMAN & GREEN JPCA 2002, 106, 11141-11149",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCtCds",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Ct 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        """,
    node="Cds-Cds(Cds-Cds)Ct",
    index=226,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCt(Cds-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Ct 0 {1,S}
        4     CO 0 {1,S}
        """,
    node="Cds-Cds(Cds-Od)(Cds-Cds)",
    index=227,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCt(Cds-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Ct 0 {1,S}
        4     Cd 0 {1,S}
        """,
    node="Cds-Cds(Cds-Cds)Ct",
    index=228,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cds)Ct",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Ct 0 {1,S}
        5     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.89,5.26,5.98,6.37,6.67,6.78,6.89],"cal/(mol*K)"),
        H298=(7.54,"kcal/mol"),
        S298=(-14.65,"cal/(mol*K)"),
    ),
    index=229,
    short_comment="Cd-CtCd RAMAN & GREEN JPCA 2002, 106, 11141-11149",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cdd)Ct",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Ct 0 {1,S}
        5     Cdd 0 {3,D}
        """,
    node="Cds-Cds(Cds-Cdd-Cd)Ct",
    index=230,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cdd-Od)Ct",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Ct 0 {1,S}
        5     Cdd 0 {3,D} {6,D}
        6     Od 0 {5,D}
        """,
    node="Cds-Cds(Cds-Cds)(Cds-Cdd-Od)",
    index=231,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cdd-Cd)Ct",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Ct 0 {1,S}
        5     Cdd 0 {3,D} {6,D}
        6     C 0 {5,D}
        """,
    node="Cds-Cds(Cds-Cds)Ct",
    index=232,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCtCt",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.23,4.59,5.41,5.93,6.48,6.74,7.02],"cal/(mol*K)"),
        H298=(8.81,"kcal/mol"),
        S298=(-13.51,"cal/(mol*K)"),
    ),
    index=233,
    short_comment="Cd-CtCt RAMAN & GREEN JPCA 2002, 106, 11141-11149",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCbCs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.4,5.37,5.93,6.18,6.5,6.62,6.72],"cal/(mol*K)"),
        H298=(8.64,"kcal/mol"),
        S298=(-14.6,"cal/(mol*K)"),
    ),
    index=234,
    short_comment="Cd-CbCs BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCbCds",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cb 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        """,
    node="Cds-Cds(Cds-Cds)Cb",
    index=235,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCb(Cds-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cb 0 {1,S}
        4     CO 0 {1,S}
        """,
    node="Cds-Cds(Cds-Od)(Cds-Cds)",
    index=236,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cd)Cb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S}
        4     Cb 0 {1,S}
        """,
    node="Cds-Cds(Cds-Cds)Cb",
    index=237,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cds)Cb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cb 0 {1,S}
        5     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.7,6.13,6.87,7.1,7.2,7.16,7.06],"cal/(mol*K)"),
        H298=(7.18,"kcal/mol"),
        S298=(-16.5,"cal/(mol*K)"),
    ),
    index=238,
    short_comment="Cd-CbCd BOZZELLI =3D Cd/Cs/Cb + (Cd/Cs/Cd - Cd/Cs/Cs)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cdd)Cb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cb 0 {1,S}
        5     Cdd 0 {3,D}
        """,
    node="Cds-Cds(Cds-Cdd-Cd)Cb",
    index=239,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cdd-Od)Cb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cb 0 {1,S}
        5     Cdd 0 {3,D} {6,D}
        6     Od 0 {5,D}
        """,
    node="Cds-Cds(Cds-Cds)(Cds-Cdd-Od)",
    index=240,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-Cds(Cds-Cdd-Cd)Cb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cd 0 {1,S} {5,D}
        4     Cb 0 {1,S}
        5     Cdd 0 {3,D} {6,D}
        6     C 0 {5,D}
        """,
    node="Cds-Cds(Cds-Cds)Cb",
    index=241,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCbCt",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([2.22,3.14,4.54,4.11,5.06,5.79,6.71],"cal/(mol*K)"),
        H298=(6.7,"kcal/mol"),
        S298=(-17.04,"cal/(mol*K)"),
    ),
    index=242,
    short_comment="Cd-CbCt Hf=3D est S,Cp mopac nov99",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CdsCbCb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cd 0 {1,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.7,6.13,6.87,7.1,7.2,7.16,7.06],"cal/(mol*K)"),
        H298=(8,"kcal/mol"),
        S298=(-16.5,"cal/(mol*K)"),
    ),
    index=243,
    short_comment="Cd-CbCb BOZZELLI =3D Cd/Cs/Cb + (Cd/Cs/Cb - Cd/Cs/Cs)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddCsCs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)CsCs",
    index=244,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)CsCs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Od 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([9.82,10.74,11.53,13.22,13.46,14.28,15.35],"cal/(mol*K)"),
        H298=(-1.644,"kcal/mol"),
        S298=(20.02,"cal/(mol*K)"),
    ),
    index=245,
    short_comment="{CCO/C2} RAMAN & GREEN JPCA 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)CsCs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-CdsCsCs",
    index=246,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddCdsCs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     {Cd,CO} 0 {1,S}
        4     Cs 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cds)Cs",
    index=247,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Od)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     CO 0 {1,S}
        4     Cs 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)CsCs",
    index=248,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cd)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S}
        4     Cs 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)Cs",
    index=249,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cds)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     Od 0 {2,D}
        6     Cd 0 {3,D}
        """,
    node="Cds-(Cdd-Od)CsCs",
    index=250,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cdd-Cd)Cs",
    index=251,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd-Od)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([10.1,11.24,12.12,12.84,14,14.75,15.72],"cal/(mol*K)"),
        H298=(-2.07,"kcal/mol"),
        S298=(19.65,"cal/(mol*K)"),
    ),
    index=252,
    short_comment="{CCO/C/CCO} RAMAN & GREEN JPCA 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd-Cd)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)Cs",
    index=253,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cd)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S}
        4     Cs 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cds)Cs",
    index=254,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cds)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     C 0 {2,D}
        6     Cd 0 {3,D}
        """,
    node="Cds-Cds(Cds-Cds)Cs",
    index=255,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     C 0 {2,D}
        6     Cdd 0 {3,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cdd-Cd)Cs",
    index=256,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd-Od)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     C 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cds-Cds(Cds-Cdd-Od)Cs",
    index=257,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd-Cd)Cs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     C 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cds)Cs",
    index=258,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddCdsCds",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     {Cd,CO} 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cds)(Cds-Cds)",
    index=259,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Od)(Cds-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     CO 0 {1,S}
        4     CO 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)CsCs",
    index=260,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cd)(Cds-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S}
        4     CO 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)(Cds-Od)",
    index=261,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cds)(Cds-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     CO 0 {1,S}
        5     Od 0 {2,D}
        6     Cd 0 {3,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Od)Cs",
    index=262,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd)(Cds-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     CO 0 {1,S}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cdd-Cd)(Cds-Od)",
    index=263,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd-Od)(Cds-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     CO 0 {1,S}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cdd-Od)Cs",
    index=264,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd-Cd)(Cds-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     CO 0 {1,S}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)(Cds-Od)",
    index=265,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cd)(Cds-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S}
        4     Cd 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)(Cds-Cds)",
    index=266,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cds)(Cds-Cds)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Od 0 {2,D}
        6     Cd 0 {3,D}
        7     Cd 0 {4,D}
        """,
    node="Cds-(Cdd-Od)CsCs",
    index=267,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd)(Cds-Cds)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D}
        7     Cd 0 {4,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cdd-Cd)(Cds-Cds)",
    index=268,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd-Od)(Cds-Cds)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D} {8,D}
        7     Cd 0 {4,D}
        8     Od 0 {6,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cdd-Od)Cs",
    index=269,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd-Cd)(Cds-Cds)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D} {8,D}
        7     Cd 0 {4,D}
        8     C 0 {6,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)(Cds-Cds)",
    index=270,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd)(Cds-Cdd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D}
        7     Cdd 0 {4,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    index=271,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)(Cds-Cds)",
    index=272,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cdd-Od)(Cds-Cds)",
    index=273,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)(Cds-Cds)",
    index=274,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Od)(Cds-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     CO 0 {1,S}
        4     CO 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-Cds(Cds-Od)(Cds-Od)",
    index=275,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Od)(Cds-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     CO 0 {1,S}
        4     Cd 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Od)(Cds-Cds)",
    index=276,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Od)(Cds-Cds)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     C 0 {2,D}
        6     Cd 0 {4,D}
        """,
    node="Cds-Cds(Cds-Od)(Cds-Cds)",
    index=277,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Od)(Cds-Cdd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     C 0 {2,D}
        6     Cdd 0 {4,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Od)(Cds-Cdd-Cd)",
    index=278,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Od)(Cds-Cdd-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     C 0 {2,D}
        6     Cdd 0 {4,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cds-Cds(Cds-Od)(Cds-Cdd-Od)",
    index=279,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Od)(Cds-Cdd-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     C 0 {2,D}
        6     Cdd 0 {4,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Od)(Cds-Cds)",
    index=280,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cd)(Cds-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S}
        4     Cd 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cds)(Cds-Cds)",
    index=281,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cds)(Cds-Cds)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     C 0 {2,D}
        6     Cd 0 {3,D}
        7     Cd 0 {4,D}
        """,
    node="Cds-Cds(Cds-Cds)(Cds-Cds)",
    index=282,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd)(Cds-Cds)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     C 0 {2,D}
        6     Cdd 0 {3,D}
        7     Cd 0 {4,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cds)",
    index=283,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd-Od)(Cds-Cds)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     C 0 {2,D}
        6     Cdd 0 {3,D} {8,D}
        7     Cd 0 {4,D}
        8     Od 0 {6,D}
        """,
    node="Cds-Cds(Cds-Cds)(Cds-Cdd-Od)",
    index=284,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cds)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     C 0 {2,D}
        6     Cdd 0 {3,D} {8,D}
        7     Cd 0 {4,D}
        8     C 0 {6,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cds)(Cds-Cds)",
    index=285,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd)(Cds-Cdd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     C 0 {2,D}
        6     Cdd 0 {3,D}
        7     Cdd 0 {4,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    index=286,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     C 0 {2,D}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cds-Cds(Cds-Cdd-Od)(Cds-Cdd-Od)",
    index=287,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd-Od)(Cds-Cdd-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     C 0 {2,D}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cdd-Od)(Cds-Cds)",
    index=288,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     C 0 {2,D}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cds)(Cds-Cds)",
    index=289,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddCtCs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)CtCs",
    index=290,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)CtCs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)Cs",
    index=291,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)CtCs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-CdsCtCs",
    index=292,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddCtCds",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     Ct 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cds)Ct",
    index=293,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Od)Ct",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     CO 0 {1,S}
        4     Ct 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)(Cds-Od)",
    index=294,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cd)Ct",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S}
        4     Ct 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)Ct",
    index=295,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cds)Ct",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     Od 0 {2,D}
        6     Cd 0 {3,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)(Cds-Cds)",
    index=296,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd)Ct",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cdd-Cd)Ct",
    index=297,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd-Od)Ct",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cdd-Od)(Cds-Cds)",
    index=298,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd-Cd)Ct",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)Ct",
    index=299,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cd)Ct",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S}
        4     Ct 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cds)Ct",
    index=300,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cds)Ct",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     C 0 {2,D}
        6     Cd 0 {3,D}
        """,
    node="Cds-Cds(Cds-Cds)Ct",
    index=301,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd)Ct",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     C 0 {2,D}
        6     Cdd 0 {3,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cdd-Cd)Ct",
    index=302,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd-Od)Ct",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     C 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cds-Cds(Cds-Cdd-Od)Ct",
    index=303,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd-Cd)Ct",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     C 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cds)Ct",
    index=304,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddCtCt",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)CtCt",
    index=305,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)CtCt",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)(Cds-Cds)",
    index=306,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)CtCt",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-CdsCtCt",
    index=307,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddCbCs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)CbCs",
    index=308,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)CbCs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)Cs",
    index=309,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)CbCs",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-CdsCbCs",
    index=310,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddCbCds",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     Cb 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cds)Cb",
    index=311,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Od)Cb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     CO 0 {1,S}
        4     Cb 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)(Cds-Od)",
    index=312,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cd)Cb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S}
        4     Cb 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)Cb",
    index=313,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cds)Cb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Od 0 {2,D}
        6     Cd 0 {3,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)(Cds-Cds)",
    index=314,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd)Cb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cdd-Cd)Cb",
    index=315,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd-Od)Cb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cdd-Od)(Cds-Cds)",
    index=316,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)(Cds-Cdd-Cd)Cb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Od 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)Cb",
    index=317,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cd)Cb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S}
        4     Cb 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cds)Cb",
    index=318,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cds)Cb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     C 0 {2,D}
        6     Cd 0 {3,D}
        """,
    node="Cds-Cds(Cds-Cds)Cb",
    index=319,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd)Cb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     C 0 {2,D}
        6     Cdd 0 {3,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cdd-Cd)Cb",
    index=320,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd-Od)Cb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     C 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cds-Cds(Cds-Cdd-Od)Cb",
    index=321,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)(Cds-Cdd-Cd)Cb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     C 0 {2,D}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cds-(Cdd-Cd)(Cds-Cds)Cb",
    index=322,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddCbCt",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)CbCt",
    index=323,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)CbCt",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)Ct",
    index=324,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)CbCt",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-CdsCbCt",
    index=325,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-CddCbCb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        """,
    node="Cds-(Cdd-Cd)CbCb",
    index=326,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Od)CbCb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Od 0 {2,D}
        """,
    node="Cds-(Cdd-Od)(Cds-Cds)(Cds-Cds)",
    index=327,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cds-(Cdd-Cd)CbCb",
    group=
        """
        1  *  C 0 {2,D} {3,S} {4,S}
        2     Cdd 0 {1,D} {5,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     C 0 {2,D}
        """,
    node="Cds-CdsCbCb",
    index=328,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs",
    group=
        """
        1  *  Cs 0
        """,
    node="Cs-CsCsCsCs",
    index=329,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-HHHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     H 0 {1,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([8.43,9.84,11.14,12.41,15,17.25,20.63],"cal/(mol*K)"),
        H298=(-17.9,"kcal/mol"),
        S298=(49.41,"cal/(mol*K)"),
    ),
    index=330,
    short_comment="CHEMKIN DATABASE S(group) = S(CH4) + Rln(12)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CHHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     C 0 {1,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-CsHHH",
    index=331,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CsHHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cs 0 {1,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([6.19,7.84,9.4,10.79,13.02,14.77,17.58],"cal/(mol*K)"),
        H298=(-10.2,"kcal/mol"),
        S298=(30.41,"cal/(mol*K)"),
    ),
    index=332,
    short_comment="Cs-CsHHH BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsHHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     {Cd,CO} 0 {1,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)HHH",
    index=333,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)HHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([6.19,7.84,9.4,10.79,13.02,14.77,17.58],"cal/(mol*K)"),
        H298=(-10.08,"kcal/mol"),
        S298=(30.41,"cal/(mol*K)"),
    ),
    index=334,
    short_comment="Cs-COHHH BENSON: Cp1500 =3D Cp1000*(Cp1500/Cp1000: C/Cd/H3)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)HHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)HHH",
    index=335,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)HHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([6.19,7.84,9.4,10.79,13.02,14.77,17.58],"cal/(mol*K)"),
        H298=(-10.2,"kcal/mol"),
        S298=(30.41,"cal/(mol*K)"),
    ),
    index=336,
    short_comment="Cs-CdHHH BENSON (Assigned Cs-CsHHH)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)HHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)HHH",
    index=337,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)HHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([6.19,7.84,9.4,10.79,13.02,14.77,17.58],"cal/(mol*K)"),
        H298=(-10.08,"kcal/mol"),
        S298=(30.41,"cal/(mol*K)"),
    ),
    index=338,
    short_comment="{CCO/C/H3} RAMAN & GREEN JPCA 2002, 106, 7937-7949, assigened same value as Cs-CsHHH",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)HHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)HHH",
    index=339,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtHHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Ct 0 {1,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([6.19,7.84,9.4,10.79,13.02,14.77,17.58],"cal/(mol*K)"),
        H298=(-10.2,"kcal/mol"),
        S298=(30.41,"cal/(mol*K)"),
    ),
    index=340,
    short_comment="Cs-CtHHH BENSON (Assigned Cs-CsHHH)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbHHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([6.19,7.84,9.4,10.79,13.02,14.77,17.58],"cal/(mol*K)"),
        H298=(-10.2,"kcal/mol"),
        S298=(30.41,"cal/(mol*K)"),
    ),
    index=341,
    short_comment="Cs-CbHHH BENSON (Assigned Cs-CsHHH)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-OsHHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Os 0 {1,S}
        3     H 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([6.19,7.84,9.4,10.79,13.03,14.77,17.58],"cal/(mol*K)"),
        H298=(-10.1,"kcal/mol"),
        S298=(30.41,"cal/(mol*K)"),
    ),
    index=342,
    short_comment="Cs-OHHH BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-OsOsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Os 0 {1,S}
        3     Os 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.5,6.95,8.25,9.35,11.07,12.34,12.34],"cal/(mol*K)"),
        H298=(-15.23,"kcal/mol"),
        S298=(9.42,"cal/(mol*K)"),
    ),
    index=343,
    short_comment="Cs-OOHH PEDLEY Hf, BOZZELLI C/C2/H2 !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-OsOsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Os 0 {1,S}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.54,6,7.17,8.05,9.31,10.05,10.05],"cal/(mol*K)"),
        H298=(-21.23,"kcal/mol"),
        S298=(-12.07,"cal/(mol*K)"),
    ),
    index=344,
    short_comment="Cs-OOOH BOZZELLI del C/C2/O - C/C3/O, series !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CCHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     C 0 {1,S}
        3     C 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-CsCsHH",
    index=345,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CsCsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cs 0 {1,S}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.5,6.95,8.25,9.35,11.07,12.34,14.25],"cal/(mol*K)"),
        H298=(-4.93,"kcal/mol"),
        S298=(9.42,"cal/(mol*K)"),
    ),
    index=346,
    short_comment="Cs-CsCsHH BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsCsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     {Cd,CO} 0 {1,S}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CsHH",
    index=347,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([6.2,7.7,8.7,9.5,11.1,12.2,14.07],"cal/(mol*K)"),
        H298=(-5.2,"kcal/mol"),
        S298=(9.6,"cal/(mol*K)"),
    ),
    index=348,
    short_comment="Cs-COCsHH BENSON Cp1500 =3D Cp1000*(Cp1500/Cp1000: C/C/Cd/H2)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CsHH",
    index=349,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.12,6.86,8.32,9.49,11.22,12.48,14.36],"cal/(mol*K)"),
        H298=(-4.76,"kcal/mol"),
        S298=(9.8,"cal/(mol*K)"),
    ),
    index=350,
    short_comment="Cs-CdCsHH BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CsHH",
    index=351,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.35,6.83,8.25,9.45,11.19,12.46,14.34],"cal/(mol*K)"),
        H298=(-5.723,"kcal/mol"),
        S298=(9.37,"cal/(mol*K)"),
    ),
    index=352,
    short_comment="{C/C/H2/CCO} RAMAN & GREEN JPCA 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CsHH",
    index=353,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsCdsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     {Cd,CO} 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)HH",
    index=354,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)HH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.03,7.44,9.16,10.49,12.17,13.57,13.57],"cal/(mol*K)"),
        H298=(-7.6,"kcal/mol"),
        S298=(5.82,"cal/(mol*K)"),
    ),
    index=355,
    short_comment="Cs-COCOHH BENSON Hf, Mopac =3D S,Cp nov99 !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)HH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)HH",
    index=356,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)HH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.75,7.11,8.92,10.32,12.16,13.61,13.61],"cal/(mol*K)"),
        H298=(-3.8,"kcal/mol"),
        S298=(6.31,"cal/(mol*K)"),
    ),
    index=357,
    short_comment="Cs-COCdHH BENSON Hf, Mopac =3D S,Cp nov99 !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)HH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)HH",
    index=358,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)HH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsHH",
    index=359,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)HH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)HH",
    index=360,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)HH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)HH",
    index=361,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)HH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.7,6.8,8.4,9.6,11.3,12.6,14.4],"cal/(mol*K)"),
        H298=(-4.29,"kcal/mol"),
        S298=(10.2,"cal/(mol*K)"),
    ),
    index=362,
    short_comment="Cs-CdCdHH BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cds)HH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cds)HH",
    index=363,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cds)HH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsHH",
    index=364,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cds)HH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)HH",
    index=365,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)HH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)HH",
    index=366,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)HH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([6.68,8.28,9.58,10.61,12.04,13.13,14.87],"cal/(mol*K)"),
        H298=(-5.301,"kcal/mol"),
        S298=(7.18,"cal/(mol*K)"),
    ),
    index=367,
    short_comment="{C/H2/CCO2} RAMAN & GREEN JPCA 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)HH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)HH",
    index=368,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)HH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)HH",
    index=369,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Ct 0 {1,S}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.95,6.56,7.93,9.08,10.86,12.19,14.2],"cal/(mol*K)"),
        H298=(-4.73,"kcal/mol"),
        S298=(10.3,"cal/(mol*K)"),
    ),
    index=370,
    short_comment="Cs-CtCsHH BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCdsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Ct 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtHH",
    index=371,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CtHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Ct 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.85,6.22,8.01,9.43,11.29,12.76,12.76],"cal/(mol*K)"),
        H298=(-5.4,"kcal/mol"),
        S298=(7.68,"cal/(mol*K)"),
    ),
    index=372,
    short_comment="Cs-COCtHH BENSON Hf, Mopac =3D S,Cp nov99 !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CtHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Ct 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtHH",
    index=373,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CtHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.4,6.33,7.9,9.16,10.93,12.29,13.43],"cal/(mol*K)"),
        H298=(-3.49,"kcal/mol"),
        S298=(9.31,"cal/(mol*K)"),
    ),
    index=374,
    short_comment="Cs-CtCdHH RAMAN & GREEN JPCA 2002, 106, 11141-11149",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CtHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CtHH",
    index=375,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CtHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)HH",
    index=376,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CtHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CtHH",
    index=377,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCtHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4,6.07,7.71,9.03,10.88,12.3,12.48],"cal/(mol*K)"),
        H298=(-0.82,"kcal/mol"),
        S298=(10.04,"cal/(mol*K)"),
    ),
    index=378,
    short_comment="Cs-CtCtHH RAMAN & GREEN JPCA 2002, 106, 11141-11149",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Cs 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.84,7.61,8.98,10.01,11.49,12.54,13.76],"cal/(mol*K)"),
        H298=(-4.86,"kcal/mol"),
        S298=(9.34,"cal/(mol*K)"),
    ),
    index=379,
    short_comment="Cs-CbCsHH BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCdsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbHH",
    index=380,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CbHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cb 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.38,7.59,9.25,10.51,12.19,13.52,13.52],"cal/(mol*K)"),
        H298=(-5.4,"kcal/mol"),
        S298=(5.89,"cal/(mol*K)"),
    ),
    index=381,
    short_comment="Cs-COCbHH BENSON Hf, Mopac =3D S,Cp nov99 !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CbHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cb 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbHH",
    index=382,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CbHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.51,6.76,8.61,10.01,11.97,13.4,15.47],"cal/(mol*K)"),
        H298=(-4.29,"kcal/mol"),
        S298=(2,"cal/(mol*K)"),
    ),
    index=383,
    short_comment="Cs-CbCdHH Hf=Stein S,Cp=3D mopac nov99",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CbHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CbHH",
    index=384,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CbHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)HH",
    index=385,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CbHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CbHH",
    index=386,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCtHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.28,6.43,8.16,9.5,11.36,12.74,13.7],"cal/(mol*K)"),
        H298=(-4.29,"kcal/mol"),
        S298=(9.84,"cal/(mol*K)"),
    ),
    index=387,
    short_comment="Cs-CbCtHH Hf=Stein S,Cp=3D mopac nov99",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.67,7.7,9.31,10.52,12.21,13.47,15.11],"cal/(mol*K)"),
        H298=(-4.29,"kcal/mol"),
        S298=(8.07,"cal/(mol*K)"),
    ),
    index=388,
    short_comment="Cs-CbCbHH Hf=3Dbsn/Cs/Cd2/H2 S,Cp=3D mopac nov99",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CCCH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     C 0 {1,S}
        3     C 0 {1,S}
        4     C 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-CsCsCsH",
    index=389,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CsCsCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cs 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.54,6,7.17,8.05,9.31,10.05,11.17],"cal/(mol*K)"),
        H298=(-1.9,"kcal/mol"),
        S298=(-12.07,"cal/(mol*K)"),
    ),
    index=390,
    short_comment="Cs-CsCsCsH BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsCsCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     {Cd,CO} 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CsCsH",
    index=391,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CsCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.16,5.91,7.34,8.19,9.46,10.19,10.19],"cal/(mol*K)"),
        H298=(-1.7,"kcal/mol"),
        S298=(-11.7,"cal/(mol*K)"),
    ),
    index=392,
    short_comment="Cs-COCsCsH BOZZELLI - BENSON Hf, S, -C/C2Cd/H =3D Cp !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CsCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CsCsH",
    index=393,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CsCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.16,5.91,7.34,8.19,9.46,10.19,11.28],"cal/(mol*K)"),
        H298=(-1.48,"kcal/mol"),
        S298=(-11.69,"cal/(mol*K)"),
    ),
    index=394,
    short_comment="Cs-CdCsCsH BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CsCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CsCsH",
    index=395,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CsCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.96,6.35,7.61,8.54,9.65,10.35,11.19],"cal/(mol*K)"),
        H298=(-3.634,"kcal/mol"),
        S298=(-12.31,"cal/(mol*K)"),
    ),
    index=396,
    short_comment="{C/C2/CCO/H} RAMAN & GREEN JPCA 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CsCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.16,5.91,7.34,8.19,9.46,10.19,11.28],"cal/(mol*K)"),
        H298=(-1.48,"kcal/mol"),
        S298=(-11.69,"cal/(mol*K)"),
    ),
    index=397,
    short_comment="Cs-CdCsCsH BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCsCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Ct 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.99,5.61,6.85,7.78,9.1,9.9,11.12],"cal/(mol*K)"),
        H298=(-1.72,"kcal/mol"),
        S298=(-11.19,"cal/(mol*K)"),
    ),
    index=398,
    short_comment="Cs-CtCsCsH BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCsCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.88,6.66,7.9,8.75,9.73,10.25,10.68],"cal/(mol*K)"),
        H298=(-0.98,"kcal/mol"),
        S298=(-12.15,"cal/(mol*K)"),
    ),
    index=399,
    short_comment="Cs-CbCsCsH BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsCdsCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     {Cd,CO} 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CsH",
    index=400,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)CsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-CsCsCsH",
    index=401,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)CsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CsH",
    index=402,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)CsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Od)CsCsH",
    index=403,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)CsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)CsH",
    index=404,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)CsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsCsH",
    index=405,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)CsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CsH",
    index=406,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)CsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CsH",
    index=407,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)CsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.28,6.54,7.67,8.48,9.45,10.18,11.24],"cal/(mol*K)"),
        H298=(-1.1,"kcal/mol"),
        S298=(-13.03,"cal/(mol*K)"),
    ),
    index=408,
    short_comment="Cs-CdCdCsH RAMAN & GREEN JPCA 2002, 106, 11141-11149",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cds)CsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cds)CsH",
    index=409,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cds)CsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsCsH",
    index=410,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cds)CsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CsH",
    index=411,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)CsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CsH",
    index=412,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([6.33,7.96,9.13,9.91,10.7,11.19,11.81],"cal/(mol*K)"),
        H298=(-3.714,"kcal/mol"),
        S298=(-14.12,"cal/(mol*K)"),
    ),
    index=413,
    short_comment="{C/C/H/CCO2} RAMAN & GREEN JPCA 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CsH",
    index=414,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CsH",
    index=415,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCdsCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtCsH",
    index=416,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CtCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CsH",
    index=417,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CtCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtCsH",
    index=418,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CtCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.55,7.21,8.39,9.17,10,10.61,10.51],"cal/(mol*K)"),
        H298=(-6.9,"kcal/mol"),
        S298=(-13.48,"cal/(mol*K)"),
    ),
    index=419,
    short_comment="Cs-CtCdCsH RAMAN & GREEN JPCA 2002, 106, 11141-11149",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CtCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CtCsH",
    index=420,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CtCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CsH",
    index=421,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CtCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CtCsH",
    index=422,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCdsCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCsH",
    index=423,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CbCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CsH",
    index=424,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CbCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCsH",
    index=425,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CbCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.5,6.57,8.07,8.89,9.88,10.39,10.79],"cal/(mol*K)"),
        H298=(-1.56,"kcal/mol"),
        S298=(-11.77,"cal/(mol*K)"),
    ),
    index=426,
    short_comment="Cs-CbCdCsH BOZZELLI =3D Cs/Cs2/Cd/H + (Cs/Cs2/Cb/H - Cs/Cs3/H)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CbCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CbCsH",
    index=427,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CbCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CsH",
    index=428,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CbCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CbCsH",
    index=429,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCtCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.27,5.32,6.9,8.03,9.33,10.21,9.38],"cal/(mol*K)"),
        H298=(1.72,"kcal/mol"),
        S298=(-11.61,"cal/(mol*K)"),
    ),
    index=430,
    short_comment="Cs-CtCtCsH RAMAN & GREEN JPCA 2002, 106, 11141-11149",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCtCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.33,6.27,7.58,8.48,9.52,10.1,10.63],"cal/(mol*K)"),
        H298=(-1.55,"kcal/mol"),
        S298=(-11.65,"cal/(mol*K)"),
    ),
    index=431,
    short_comment="Cs-CbCtCsH BOZZELLI =3D Cs/Cs2/Cb/H + (Cs/Cs2/Ct/H - Cs/Cs3/H)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbCsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.22,7.32,8.63,8.45,10.15,10.45,10.89],"cal/(mol*K)"),
        H298=(-1.06,"kcal/mol"),
        S298=(-12.23,"cal/(mol*K)"),
    ),
    index=432,
    short_comment="Cs-CbCbCsCs BOZZELLI =3D Cs/Cs2/Cb/H + (Cs/Cs2/Cb/H - Cs/Cs3/H)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsCdsCdsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     {Cd,CO} 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H",
    index=433,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Od)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     CO 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-CsCsCsH",
    index=434,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)H",
    index=435,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     H 0 {1,S}
        6     Cd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)CsH",
    index=436,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     H 0 {1,S}
        6     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)H",
    index=437,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     H 0 {1,S}
        6     Cdd 0 {4,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsCsH",
    index=438,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     H 0 {1,S}
        6     Cdd 0 {4,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)H",
    index=439,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)(Cds-Cd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Cd 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)H",
    index=440,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     H 0 {1,S}
        6     Cd 0 {3,D}
        7     Cd 0 {4,D}
        """,
    node="Cs-(Cds-Od)CsCsH",
    index=441,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)(Cds-Cds)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     H 0 {1,S}
        6     Cdd 0 {3,D}
        7     Cd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cds)H",
    index=442,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     H 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cd 0 {4,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Od)CsH",
    index=443,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cds)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     H 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cd 0 {4,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)H",
    index=444,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)(Cds-Cdd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     H 0 {1,S}
        6     Cdd 0 {3,D}
        7     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)H",
    index=445,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     H 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CsH",
    index=446,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     H 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)H",
    index=447,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     H 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)H",
    index=448,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)(Cds-Cd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Cd 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H",
    index=449,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cd 0 {4,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.51,5.96,7.13,7.98,9.06,9.9,11.23],"cal/(mol*K)"),
        H298=(0.41,"kcal/mol"),
        S298=(-11.82,"cal/(mol*K)"),
    ),
    index=450,
    short_comment="Cs-CdCdCdH RAMAN & GREEN JPC 2002",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)H",
    index=451,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cdd 0 {4,D} {9,D}
        9     Od 0 {8,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsCsH",
    index=452,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cdd 0 {4,D} {9,D}
        9     C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H",
    index=453,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd)(Cds-Cdd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D}
        8     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)H",
    index=454,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D} {9,D}
        8     Cdd 0 {4,D} {10,D}
        9     Od 0 {7,D}
        10    Od 0 {8,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CsH",
    index=455,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Cd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D} {9,D}
        8     Cdd 0 {4,D} {10,D}
        9     Od 0 {7,D}
        10    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)H",
    index=456,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D} {9,D}
        8     Cdd 0 {4,D} {10,D}
        9     C 0 {7,D}
        10    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H",
    index=457,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        8     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)H",
    index=458,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     Od 0 {6,D}
        10    Od 0 {7,D}
        11    Od 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H",
    index=459,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     Od 0 {6,D}
        10    Od 0 {7,D}
        11    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)H",
    index=460,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     Od 0 {6,D}
        10    C 0 {7,D}
        11    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)H",
    index=461,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     C 0 {6,D}
        10    C 0 {7,D}
        11    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H",
    index=462,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCdsCdsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtH",
    index=463,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)CtH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Ct 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)H",
    index=464,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)CtH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Ct 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CtH",
    index=465,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)CtH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)H",
    index=466,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)CtH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)CtH",
    index=467,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)CtH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)H",
    index=468,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)CtH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CtH",
    index=469,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)CtH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Ct 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtH",
    index=470,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)CtH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([6.68,7.85,8.62,9.16,9.81,10.42,10.49],"cal/(mol*K)"),
        H298=(1.88,"kcal/mol"),
        S298=(-13.75,"cal/(mol*K)"),
    ),
    index=471,
    short_comment="Cs-CtCdCdH RAMAN & GREEN JPCA 2002, 106, 11141-11149",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cds)CtH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cds)CtH",
    index=472,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cds)CtH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)H",
    index=473,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cds)CtH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtH",
    index=474,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)CtH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtH",
    index=475,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CtH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)H",
    index=476,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CtH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CtH",
    index=477,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtH",
    index=478,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCdsCdsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbH",
    index=479,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)CbH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cb 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)H",
    index=480,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)CbH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Cb 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CbH",
    index=481,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)CbH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)H",
    index=482,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)CbH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)CbH",
    index=483,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)CbH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)H",
    index=484,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)CbH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CbH",
    index=485,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)CbH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Cb 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbH",
    index=486,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)CbH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.12,6.51,8.24,9,10.03,10.53,10.89],"cal/(mol*K)"),
        H298=(-1.39,"kcal/mol"),
        S298=(-11.39,"cal/(mol*K)"),
    ),
    index=487,
    short_comment="Cs-CbCdCdH BOZZELLI =3D Cs/Cs/Cd2/H + (Cs/Cs2/Cb/H - Cs/Cs3/H)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cds)CbH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cds)CbH",
    index=488,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cds)CbH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)H",
    index=489,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cds)CbH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbH",
    index=490,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)CbH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbH",
    index=491,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CbH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)H",
    index=492,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CbH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CbH",
    index=493,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbH",
    index=494,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCtCdsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-CtCt(Cds-Cds)H",
    index=495,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCt(Cds-Od)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     CO 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)H",
    index=496,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCt(Cds-Cd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     Cd 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-CtCt(Cds-Cds)H",
    index=497,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCt(Cds-Cds)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     H 0 {1,S}
        6     Cd 0 {4,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.58,5.68,7.11,8.12,9.27,10.13,9.44],"cal/(mol*K)"),
        H298=(4.73,"kcal/mol"),
        S298=(-11.46,"cal/(mol*K)"),
    ),
    index=498,
    short_comment="Cs-CtCtCdH RAMAN & GREEN JPCA 2002, 106, 11141-11149",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCt(Cds-Cdd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     H 0 {1,S}
        6     Cdd 0 {4,D}
        """,
    node="Cs-CtCt(Cds-Cdd-Cd)H",
    index=499,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCt(Cds-Cdd-Od)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     H 0 {1,S}
        6     Cdd 0 {4,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)H",
    index=500,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCt(Cds-Cdd-Cd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     H 0 {1,S}
        6     Cdd 0 {4,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-CtCt(Cds-Cds)H",
    index=501,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCtCdsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-CbCt(Cds-Cds)H",
    index=502,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCt(Cds-Od)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     CO 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CtH",
    index=503,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCt(Cds-Cd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     Cd 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-CbCt(Cds-Cds)H",
    index=504,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCt(Cds-Cds)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     H 0 {1,S}
        6     Cd 0 {4,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtH",
    index=505,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCt(Cds-Cdd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     H 0 {1,S}
        6     Cdd 0 {4,D}
        """,
    node="Cs-CbCt(Cds-Cdd-Cd)H",
    index=506,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCt(Cds-Cdd-Od)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     H 0 {1,S}
        6     Cdd 0 {4,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CtH",
    index=507,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCt(Cds-Cdd-Cd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     H 0 {1,S}
        6     Cdd 0 {4,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-CbCt(Cds-Cds)H",
    index=508,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbCdsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-CbCb(Cds-Cds)H",
    index=509,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCb(Cds-Od)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     CO 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)H",
    index=510,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCb(Cds-Cd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Cd 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-CbCbCdsH",
    index=0,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCb(Cds-Cds)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     H 0 {1,S}
        6     Cd 0 {4,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)H",
    index=511,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCb(Cds-Cdd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     H 0 {1,S}
        6     Cdd 0 {4,D}
        """,
    node="Cs-CbCb(Cds-Cdd-Cd)H",
    index=512,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCb(Cds-Cdd-Od)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     H 0 {1,S}
        6     Cdd 0 {4,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)H",
    index=513,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCb(Cds-Cdd-Cd)H",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     H 0 {1,S}
        6     Cdd 0 {4,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-CbCb(Cds-Cds)H",
    index=514,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCtCtH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.03,5.27,6.78,7.88,9.14,10.08,8.47],"cal/(mol*K)"),
        H298=(10.11,"kcal/mol"),
        S298=(-10.46,"cal/(mol*K)"),
    ),
    index=515,
    short_comment="Cs-CtCtCtH RAMAN & GREEN JPCA 2002, 106, 11141-11149",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCtCtH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-CtCt(Cds-Cds)H",
    index=516,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbCtH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtH",
    index=517,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbCbH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.56,7.98,9.36,10.15,10.57,10.65,9.7],"cal/(mol*K)"),
        H298=(-0.34,"kcal/mol"),
        S298=(-12.31,"cal/(mol*K)"),
    ),
    index=518,
    short_comment="Cs-CbCbCbH BOZZELLI =3D Cs/Cs/Cb2/H + (Cs/Cs2/Cb/H - Cs/Cs3/H)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CCCC",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     C 0 {1,S}
        3     C 0 {1,S}
        4     C 0 {1,S}
        5     C 0 {1,S}
        """,
    node="Cs-CsCsCsCs",
    index=519,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CsCsCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cs 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.37,6.13,7.36,8.12,8.77,8.76,8.12],"cal/(mol*K)"),
        H298=(0.5,"kcal/mol"),
        S298=(-35.1,"cal/(mol*K)"),
    ),
    index=520,
    short_comment="Cs-CsCsCsCs BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsCsCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     {Cd,CO} 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CsCsCs",
    index=521,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CsCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.99,6.04,7.43,8.26,8.92,8.96,8.23],"cal/(mol*K)"),
        H298=(1.4,"kcal/mol"),
        S298=(-34.72,"cal/(mol*K)"),
    ),
    index=522,
    short_comment="Cs-COCsCsCs Hf BENSON S,Cp =3D C/Cd/C3",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CsCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CsCsCs",
    index=523,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CsCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cd 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.99,6.04,7.43,8.26,8.92,8.96,8.23],"cal/(mol*K)"),
        H298=(1.68,"kcal/mol"),
        S298=(-34.72,"cal/(mol*K)"),
    ),
    index=524,
    short_comment="Cs-CdCsCsCs BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CsCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CsCsCs",
    index=525,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CsCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.48,6.06,7.31,8.07,8.59,8.66,8.29],"cal/(mol*K)"),
        H298=(-2.896,"kcal/mol"),
        S298=(-34.87,"cal/(mol*K)"),
    ),
    index=526,
    short_comment="{C/C3/CCO} RAMAN & GREEN JPCA 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CsCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CsCsCs",
    index=527,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCsCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.37,6.79,8.09,8.78,9.19,8.96,7.63],"cal/(mol*K)"),
        H298=(2.81,"kcal/mol"),
        S298=(-35.18,"cal/(mol*K)"),
    ),
    index=528,
    short_comment="Cs-CtCsCsCs BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCsCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.37,6.79,8.09,8.78,9.19,8.96,7.63],"cal/(mol*K)"),
        H298=(2.81,"kcal/mol"),
        S298=(-35.18,"cal/(mol*K)"),
    ),
    index=529,
    short_comment="Cs-CbCsCsCs BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsCdsCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     {Cd,CO} 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CsCs",
    index=530,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)CsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-CsCsCsCs",
    index=531,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)CsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CsCs",
    index=532,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)CsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Od)CsCsCs",
    index=533,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)CsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)CsCs",
    index=534,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)CsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsCsCs",
    index=535,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)CsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CsCs",
    index=536,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)CsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CsCs",
    index=537,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)CsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.99,6.04,7.43,8.26,8.92,8.96,8.23],"cal/(mol*K)"),
        H298=(1.68,"kcal/mol"),
        S298=(-34.72,"cal/(mol*K)"),
    ),
    index=538,
    short_comment="Cs-CdCdCsCs BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cds)CsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cds)CsCs",
    index=539,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cds)CsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsCsCs",
    index=540,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cds)CsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CsCs",
    index=541,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)CsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CsCs",
    index=542,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([6.73,8.1,9.02,9.53,9.66,9.52,8.93],"cal/(mol*K)"),
        H298=(-2.987,"kcal/mol"),
        S298=(-36.46,"cal/(mol*K)"),
    ),
    index=543,
    short_comment="{C/C2/CCO2} RAMAN & GREEN JPCA 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CsCs",
    index=544,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CsCs",
    index=545,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCdsCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtCsCs",
    index=546,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CtCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CsCs",
    index=547,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CtCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtCsCs",
    index=548,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CtCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cd 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.99,6.7,8.16,8.92,9.34,9.16,7.14],"cal/(mol*K)"),
        H298=(2.99,"kcal/mol"),
        S298=(-34.8,"cal/(mol*K)"),
    ),
    index=549,
    short_comment="Cs-CtCdCsCs BOZZELLI =3D Cs/Cs3/Cd + (Cs/Cs3/Ct - Cs/Cs4)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CtCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CtCsCs",
    index=550,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CtCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CsCs",
    index=551,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CtCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CtCsCs",
    index=552,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCdsCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCsCs",
    index=553,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CbCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CsCs",
    index=554,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CbCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCsCs",
    index=555,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CbCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cd 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.99,6.7,8.16,8.92,9.34,9.16,7.14],"cal/(mol*K)"),
        H298=(2.99,"kcal/mol"),
        S298=(-34.8,"cal/(mol*K)"),
    ),
    index=556,
    short_comment="Cs-CbCdCsCs BOZZELLI =3D Cs/Cs3/Cb + (Cs/Cs3/Cd - Cs/Cs4)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CbCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CbCsCs",
    index=557,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CbCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CsCs",
    index=558,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CbCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CbCsCs",
    index=559,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCtCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.57,5.98,7.51,8.37,9,9.02,8.34],"cal/(mol*K)"),
        H298=(1.16,"kcal/mol"),
        S298=(-35.26,"cal/(mol*K)"),
    ),
    index=560,
    short_comment="Cs-CtCtCsCs BOZZELLI =3D Cs/Cs3/Ct + (Cs/Cs3/Ct - Cs/Cs4)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCtCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.57,5.98,7.51,8.37,9,9.02,8.34],"cal/(mol*K)"),
        H298=(1.16,"kcal/mol"),
        S298=(-35.26,"cal/(mol*K)"),
    ),
    index=561,
    short_comment="Cs-CbCtCsCs BOZZELLI =3D Cs/Cs3/Cb + (Cs/Cs3/Ct - Cs/Cs4)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.57,5.98,7.51,8.37,9,9.02,8.34],"cal/(mol*K)"),
        H298=(1.16,"kcal/mol"),
        S298=(-35.26,"cal/(mol*K)"),
    ),
    index=562,
    short_comment="Cs-CbCbCsCs BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsCdsCdsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     {Cd,CO} 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cs",
    index=563,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Od)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     CO 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-CsCsCsCs",
    index=564,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cd)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Cs",
    index=565,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Cs 0 {1,S}
        6     Cd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)CsCs",
    index=566,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Cs 0 {1,S}
        6     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)Cs",
    index=567,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Cs 0 {1,S}
        6     Cdd 0 {4,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsCsCs",
    index=568,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Cs 0 {1,S}
        6     Cdd 0 {4,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Cs",
    index=569,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)(Cds-Cd)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Cd 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Cs",
    index=570,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cs 0 {1,S}
        6     Cd 0 {3,D}
        7     Cd 0 {4,D}
        """,
    node="Cs-(Cds-Od)CsCsCs",
    index=571,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)(Cds-Cds)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cs 0 {1,S}
        6     Cdd 0 {3,D}
        7     Cd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cds)Cs",
    index=572,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cs 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cd 0 {4,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Od)CsCs",
    index=573,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cds)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cs 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cd 0 {4,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Cs",
    index=574,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)(Cds-Cdd)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cs 0 {1,S}
        6     Cdd 0 {3,D}
        7     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs",
    index=575,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cs 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CsCs",
    index=576,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cs 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)Cs",
    index=577,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cs 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Cs",
    index=578,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)(Cds-Cd)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Cd 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cs",
    index=579,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cs 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cd 0 {4,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.32,5.86,7.57,8.54,9.22,9.36,8.45],"cal/(mol*K)"),
        H298=(2.54,"kcal/mol"),
        S298=(-33.96,"cal/(mol*K)"),
    ),
    index=580,
    short_comment="Cs-CdCdCdCs BOZZELLI =3D Cs/Cs2/Cd2 + (Cs/Cs3/Cd - Cs/Cs4)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cs 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)Cs",
    index=581,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cs 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cdd 0 {4,D} {9,D}
        9     Od 0 {8,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsCsCs",
    index=582,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cs 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cdd 0 {4,D} {9,D}
        9     C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cs",
    index=583,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd)(Cds-Cdd)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cs 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D}
        8     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs",
    index=584,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cs 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D} {9,D}
        8     Cdd 0 {4,D} {10,D}
        9     Od 0 {7,D}
        10    Od 0 {8,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CsCs",
    index=585,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Cd)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cs 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D} {9,D}
        8     Cdd 0 {4,D} {10,D}
        9     Od 0 {7,D}
        10    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Cs",
    index=586,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cs 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D} {9,D}
        8     Cdd 0 {4,D} {10,D}
        9     C 0 {7,D}
        10    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cs",
    index=587,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        8     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs",
    index=588,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     Od 0 {6,D}
        10    Od 0 {7,D}
        11    Od 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)Cs",
    index=589,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     Od 0 {6,D}
        10    Od 0 {7,D}
        11    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)Cs",
    index=590,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     Od 0 {6,D}
        10    C 0 {7,D}
        11    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Cs",
    index=591,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     C 0 {6,D}
        10    C 0 {7,D}
        11    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cs",
    index=592,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCdsCdsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Ct 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtCs",
    index=593,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)CtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Cs",
    index=594,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)CtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CtCs",
    index=595,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)CtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Cs",
    index=596,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)CtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)CtCs",
    index=597,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)CtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)Cs",
    index=598,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)CtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CtCs",
    index=599,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)CtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtCs",
    index=600,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)CtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cs",
    index=601,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cds)CtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cds)CtCs",
    index=602,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cds)CtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Cs",
    index=603,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cds)CtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtCs",
    index=604,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)CtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtCs",
    index=605,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)Cs",
    index=606,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CtCs",
    index=607,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtCs",
    index=608,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCdsCdsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbCs",
    index=609,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)CbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Cs",
    index=610,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)CbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CbCs",
    index=611,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)CbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        6     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Cs",
    index=612,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)CbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)CbCs",
    index=613,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)CbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)Cs",
    index=614,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)CbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CbCs",
    index=615,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)CbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbCs",
    index=616,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)CbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cs",
    index=617,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cds)CbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cds)CbCs",
    index=618,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cds)CbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Cs",
    index=619,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cds)CbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbCs",
    index=620,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)CbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbCs",
    index=621,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)Cs",
    index=622,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CbCs",
    index=623,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbCs",
    index=624,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCtCdsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtCtCs",
    index=625,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CtCtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Cs",
    index=626,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CtCtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtCtCs",
    index=627,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CtCtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cd 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.99,7.36,8.89,9.58,9.76,9.16,7.25],"cal/(mol*K)"),
        H298=(5.1,"kcal/mol"),
        S298=(-34.88,"cal/(mol*K)"),
    ),
    index=628,
    short_comment="Cs-CtCtCdCs BOZZELLI =3D Cs/Cd2/Cs2 + (Cs/Cs3/Ct - Cs/Cs4)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CtCtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CtCtCs",
    index=629,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CtCtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Cs",
    index=630,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CtCtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CtCtCs",
    index=631,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCtCdsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCtCs",
    index=632,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CbCtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CtCs",
    index=633,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CbCtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCtCs",
    index=634,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CbCtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cd 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.99,7.36,8.89,9.58,9.76,9.16,7.25],"cal/(mol*K)"),
        H298=(5.1,"kcal/mol"),
        S298=(-34.88,"cal/(mol*K)"),
    ),
    index=635,
    short_comment="Cs-CbCtCdCs BOZZELLI =3D Cs/Cb/Cd/Cs2 + (Cs/Cs3/Ct - Cs/Cs4)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CbCtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CbCtCs",
    index=636,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CbCtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CtCs",
    index=637,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CbCtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.99,7.36,8.89,9.58,9.76,9.16,7.25],"cal/(mol*K)"),
        H298=(5.1,"kcal/mol"),
        S298=(-34.88,"cal/(mol*K)"),
    ),
    index=638,
    short_comment="Cs-CbCtCdCs BOZZELLI =3D Cs/Cb/Cd/Cs2 + (Cs/Cs3/Ct - Cs/Cs4)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbCdsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCbCs",
    index=639,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CbCbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Cs",
    index=640,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CbCbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCbCs",
    index=641,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CbCbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        6     Cd 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.99,7.36,8.89,9.58,9.76,9.16,7.25],"cal/(mol*K)"),
        H298=(5.1,"kcal/mol"),
        S298=(-34.88,"cal/(mol*K)"),
    ),
    index=642,
    short_comment="Cs-CbCbCdCs BOZZELLI =3D Cs/Cs2/Cb2 + (Cs/Cs3/Cd - Cs/Cs4)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CbCbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CbCbCs",
    index=643,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CbCbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Cs",
    index=644,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CbCbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CbCbCs",
    index=645,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCtCtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.37,8.11,9.55,10.1,10.03,9.36,6.65],"cal/(mol*K)"),
        H298=(6.23,"kcal/mol"),
        S298=(-35.34,"cal/(mol*K)"),
    ),
    index=646,
    short_comment="Cs-CtCtCtCs BOZZELLI =3D Cs/Cs2/Ct2 + (Cs/Cs3/Ct - Cs/Cs4)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCtCtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.37,8.11,9.55,10.1,10.03,9.36,6.65],"cal/(mol*K)"),
        H298=(6.23,"kcal/mol"),
        S298=(-35.34,"cal/(mol*K)"),
    ),
    index=647,
    short_comment="Cs-CbCtCtCs BOZZELLI =3D Cs/Cs2/Cb/Ct + (Cs/Cs3/Ct - Cs/Cs4)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbCtCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.37,8.11,9.55,10.1,10.03,9.36,6.65],"cal/(mol*K)"),
        H298=(6.43,"kcal/mol"),
        S298=(-35.34,"cal/(mol*K)"),
    ),
    index=648,
    short_comment="Cs-CbCbCtCs BOZZELLI =3D Cs/Cs2/Cb2 + (Cs/Cs3/Ct - Cs/Cs4)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbCbCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.37,8.11,9.55,10.1,10.03,9.36,6.65],"cal/(mol*K)"),
        H298=(6.23,"kcal/mol"),
        S298=(-35.34,"cal/(mol*K)"),
    ),
    index=649,
    short_comment="Cs-CbCbCbCs BOZZELLI =3D Cs/Cs2/Cb2 + (Cs/Cs3/Cb - Cs/Cs4)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsCdsCdsCds",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     {Cd,CO} 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     {Cd,CO} 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=650,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Od)(Cds-Od)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     CO 0 {1,S}
        5     CO 0 {1,S}
        """,
    node="Cs-CsCsCsCs",
    index=651,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Od)(Cds-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     CO 0 {1,S}
        5     Cd 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Od)(Cds-Cds)",
    index=652,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Od)(Cds-Cds)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     CO 0 {1,S}
        5     Cd 0 {1,S} {6,D}
        6     Cd 0 {5,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Od)Cs",
    index=653,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Od)(Cds-Cdd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     CO 0 {1,S}
        5     Cd 0 {1,S} {6,D}
        6     Cdd 0 {5,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)",
    index=654,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Od)(Cds-Cdd-Od)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     CO 0 {1,S}
        5     Cd 0 {1,S} {6,D}
        6     Cdd 0 {5,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsCsCs",
    index=655,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     CO 0 {1,S}
        5     Cd 0 {1,S} {6,D}
        6     Cdd 0 {5,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Od)(Cds-Cds)",
    index=656,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cd)(Cds-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S}
        5     Cd 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)(Cds-Cds)",
    index=657,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)(Cds-Cds)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Cd 0 {1,S} {7,D}
        6     Cd 0 {4,D}
        7     Cd 0 {5,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)CsCs",
    index=658,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd)(Cds-Cds)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Cd 0 {1,S} {7,D}
        6     Cdd 0 {4,D}
        7     Cd 0 {5,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)(Cds-Cds)",
    index=659,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Cd 0 {1,S} {7,D}
        6     Cdd 0 {4,D} {8,D}
        7     Cd 0 {5,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)Cs",
    index=660,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)(Cds-Cds)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Cd 0 {1,S} {7,D}
        6     Cdd 0 {4,D} {8,D}
        7     Cd 0 {5,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)(Cds-Cds)",
    index=661,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd)(Cds-Cdd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Cd 0 {1,S} {7,D}
        6     Cdd 0 {4,D}
        7     Cdd 0 {5,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    index=662,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Cd 0 {1,S} {7,D}
        6     Cdd 0 {4,D} {8,D}
        7     Cdd 0 {5,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CsCs",
    index=663,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Cd 0 {1,S} {7,D}
        6     Cdd 0 {4,D} {8,D}
        7     Cdd 0 {5,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)",
    index=664,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Cd 0 {1,S} {7,D}
        6     Cdd 0 {4,D} {8,D}
        7     Cdd 0 {5,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)(Cds-Cds)",
    index=665,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)(Cds-Cd)(Cds-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Cd 0 {1,S}
        5     Cd 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=666,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cd 0 {1,S} {8,D}
        6     Cd 0 {3,D}
        7     Cd 0 {4,D}
        8     Cd 0 {5,D}
        """,
    node="Cs-(Cds-Od)CsCsCs",
    index=667,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cdd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cd 0 {1,S} {8,D}
        6     Cd 0 {3,D}
        7     Cd 0 {4,D}
        8     Cdd 0 {5,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)",
    index=668,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cd 0 {1,S} {8,D}
        6     Cd 0 {3,D}
        7     Cd 0 {4,D}
        8     Cdd 0 {5,D} {9,D}
        9     Od 0 {8,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Od)CsCs",
    index=669,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cd 0 {1,S} {8,D}
        6     Cd 0 {3,D}
        7     Cd 0 {4,D}
        8     Cdd 0 {5,D} {9,D}
        9     C 0 {8,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=670,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)(Cds-Cdd)(Cds-Cdd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cd 0 {1,S} {8,D}
        6     Cd 0 {3,D}
        7     Cdd 0 {4,D}
        8     Cdd 0 {5,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    index=671,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cd 0 {1,S} {8,D}
        6     Cd 0 {3,D}
        7     Cdd 0 {4,D} {9,D}
        8     Cdd 0 {5,D} {10,D}
        9     Od 0 {7,D}
        10    Od 0 {8,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)Cs",
    index=672,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cd 0 {1,S} {8,D}
        6     Cd 0 {3,D}
        7     Cdd 0 {4,D} {9,D}
        8     Cdd 0 {5,D} {10,D}
        9     Od 0 {7,D}
        10    C 0 {8,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)",
    index=673,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cd 0 {1,S} {8,D}
        6     Cd 0 {3,D}
        7     Cdd 0 {4,D} {9,D}
        8     Cdd 0 {5,D} {10,D}
        9     C 0 {7,D}
        10    C 0 {8,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=674,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cd 0 {1,S} {8,D}
        6     Cdd 0 {3,D}
        7     Cdd 0 {4,D}
        8     Cdd 0 {5,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    index=675,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cd 0 {1,S} {8,D}
        6     Cdd 0 {3,D} {9,D}
        7     Cdd 0 {4,D} {10,D}
        8     Cdd 0 {5,D} {11,D}
        9     Od 0 {6,D}
        10    Od 0 {7,D}
        11    Od 0 {8,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)Cs",
    index=676,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cd 0 {1,S} {8,D}
        6     Cdd 0 {3,D} {9,D}
        7     Cdd 0 {4,D} {10,D}
        8     Cdd 0 {5,D} {11,D}
        9     Od 0 {6,D}
        10    Od 0 {7,D}
        11    C 0 {8,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    index=677,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cd 0 {1,S} {8,D}
        6     Cdd 0 {3,D} {9,D}
        7     Cdd 0 {4,D} {10,D}
        8     Cdd 0 {5,D} {11,D}
        9     Od 0 {6,D}
        10    C 0 {7,D}
        11    C 0 {8,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)",
    index=678,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cd 0 {1,S} {8,D}
        6     Cdd 0 {3,D} {9,D}
        7     Cdd 0 {4,D} {10,D}
        8     Cdd 0 {5,D} {11,D}
        9     C 0 {6,D}
        10    C 0 {7,D}
        11    C 0 {8,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=679,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)(Cds-Cd)(Cds-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Cd 0 {1,S}
        5     Cd 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=680,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cd 0 {1,S} {9,D}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cd 0 {4,D}
        9     Cd 0 {5,D}
        """,
    node="Cs-CsCsCsCs",
    index=681,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cd 0 {1,S} {9,D}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cd 0 {4,D}
        9     Cdd 0 {5,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)",
    index=682,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cd 0 {1,S} {9,D}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cd 0 {4,D}
        9     Cdd 0 {5,D} {10,D}
        10    Od 0 {9,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsCsCs",
    index=683,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cd 0 {1,S} {9,D}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cd 0 {4,D}
        9     Cdd 0 {5,D} {10,D}
        10    C 0 {9,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=684,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)(Cds-Cdd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cd 0 {1,S} {9,D}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cdd 0 {4,D}
        9     Cdd 0 {5,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    index=685,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cd 0 {1,S} {9,D}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cdd 0 {4,D} {10,D}
        9     Cdd 0 {5,D} {11,D}
        10    Od 0 {8,D}
        11    Od 0 {9,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CsCs",
    index=686,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cd 0 {1,S} {9,D}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cdd 0 {4,D} {10,D}
        9     Cdd 0 {5,D} {11,D}
        10    Od 0 {8,D}
        11    C 0 {9,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)",
    index=687,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cd 0 {1,S} {9,D}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cdd 0 {4,D} {10,D}
        9     Cdd 0 {5,D} {11,D}
        10    C 0 {8,D}
        11    C 0 {9,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=688,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cd 0 {1,S} {9,D}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D}
        8     Cdd 0 {4,D}
        9     Cdd 0 {5,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    index=689,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cd 0 {1,S} {9,D}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     Cdd 0 {5,D} {12,D}
        10    Od 0 {7,D}
        11    Od 0 {8,D}
        12    Od 0 {9,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)Cs",
    index=690,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cd 0 {1,S} {9,D}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     Cdd 0 {5,D} {12,D}
        10    Od 0 {7,D}
        11    Od 0 {8,D}
        12    C 0 {9,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    index=691,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cd 0 {1,S} {9,D}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     Cdd 0 {5,D} {12,D}
        10    Od 0 {7,D}
        11    C 0 {8,D}
        12    C 0 {9,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)",
    index=692,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cd 0 {1,S} {9,D}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     Cdd 0 {5,D} {12,D}
        10    C 0 {7,D}
        11    C 0 {8,D}
        12    C 0 {9,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=693,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cd 0 {1,S} {9,D}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        8     Cdd 0 {4,D}
        9     Cdd 0 {5,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    index=694,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cd 0 {1,S} {9,D}
        6     Cdd 0 {2,D} {10,D}
        7     Cdd 0 {3,D} {11,D}
        8     Cdd 0 {4,D} {12,D}
        9     Cdd 0 {5,D} {13,D}
        10    Od 0 {6,D}
        11    Od 0 {7,D}
        12    Od 0 {8,D}
        13    Od 0 {9,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=695,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cd 0 {1,S} {9,D}
        6     Cdd 0 {2,D} {10,D}
        7     Cdd 0 {3,D} {11,D}
        8     Cdd 0 {4,D} {12,D}
        9     Cdd 0 {5,D} {13,D}
        10    Od 0 {6,D}
        11    Od 0 {7,D}
        12    Od 0 {8,D}
        13    C 0 {9,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    index=696,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cd 0 {1,S} {9,D}
        6     Cdd 0 {2,D} {10,D}
        7     Cdd 0 {3,D} {11,D}
        8     Cdd 0 {4,D} {12,D}
        9     Cdd 0 {5,D} {13,D}
        10    Od 0 {6,D}
        11    Od 0 {7,D}
        12    C 0 {8,D}
        13    C 0 {9,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    index=697,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cd 0 {1,S} {9,D}
        6     Cdd 0 {2,D} {10,D}
        7     Cdd 0 {3,D} {11,D}
        8     Cdd 0 {4,D} {12,D}
        9     Cdd 0 {5,D} {13,D}
        10    Od 0 {6,D}
        11    C 0 {7,D}
        12    C 0 {8,D}
        13    C 0 {9,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)",
    index=698,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cd 0 {1,S} {9,D}
        6     Cdd 0 {2,D} {10,D}
        7     Cdd 0 {3,D} {11,D}
        8     Cdd 0 {4,D} {12,D}
        9     Cdd 0 {5,D} {13,D}
        10    C 0 {6,D}
        11    C 0 {7,D}
        12    C 0 {8,D}
        13    C 0 {9,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=699,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCdsCdsCds",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Ct 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     {Cd,CO} 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Ct",
    index=700,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Od)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     CO 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Od)(Cds-Cds)",
    index=701,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cd)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Ct",
    index=702,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Ct 0 {1,S}
        6     Cd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)(Cds-Cds)",
    index=703,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Ct 0 {1,S}
        6     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)Ct",
    index=704,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Ct 0 {1,S}
        6     Cdd 0 {4,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)",
    index=705,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Ct 0 {1,S}
        6     Cdd 0 {4,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Ct",
    index=706,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)(Cds-Cd)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Cd 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Ct",
    index=707,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Ct 0 {1,S}
        6     Cd 0 {3,D}
        7     Cd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=708,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)(Cds-Cds)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Ct 0 {1,S}
        6     Cdd 0 {3,D}
        7     Cd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cds)Ct",
    index=709,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Ct 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cd 0 {4,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)",
    index=710,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cds)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Ct 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cd 0 {4,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Ct",
    index=711,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)(Cds-Cdd)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Ct 0 {1,S}
        6     Cdd 0 {3,D}
        7     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct",
    index=712,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Ct 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    index=713,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Ct 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)Ct",
    index=714,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Ct 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Ct",
    index=715,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)(Cds-Cd)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Cd 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Ct",
    index=716,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Ct 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cd 0 {4,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=717,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Ct 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)Ct",
    index=718,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Ct 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cdd 0 {4,D} {9,D}
        9     Od 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)",
    index=719,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Ct 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cdd 0 {4,D} {9,D}
        9     C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Ct",
    index=720,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd)(Cds-Cdd)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Ct 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D}
        8     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct",
    index=721,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Ct 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D} {9,D}
        8     Cdd 0 {4,D} {10,D}
        9     Od 0 {7,D}
        10    Od 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    index=722,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Cd)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Ct 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D} {9,D}
        8     Cdd 0 {4,D} {10,D}
        9     Od 0 {7,D}
        10    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Ct",
    index=723,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Ct 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D} {9,D}
        8     Cdd 0 {4,D} {10,D}
        9     C 0 {7,D}
        10    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Ct",
    index=724,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        8     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct",
    index=725,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     Od 0 {6,D}
        10    Od 0 {7,D}
        11    Od 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    index=726,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     Od 0 {6,D}
        10    Od 0 {7,D}
        11    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)Ct",
    index=727,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     Od 0 {6,D}
        10    C 0 {7,D}
        11    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Ct",
    index=728,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Ct",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     C 0 {6,D}
        10    C 0 {7,D}
        11    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Ct",
    index=729,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCdsCdsCds",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     {Cd,CO} 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cb",
    index=730,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Od)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     CO 0 {1,S}
        5     Cb 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Od)(Cds-Cds)",
    index=731,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cd)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S}
        5     Cb 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Cb",
    index=732,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Cb 0 {1,S}
        6     Cd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)(Cds-Cds)",
    index=733,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Cb 0 {1,S}
        6     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)Cb",
    index=734,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Cb 0 {1,S}
        6     Cdd 0 {4,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)",
    index=735,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Cb 0 {1,S}
        6     Cdd 0 {4,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Cb",
    index=736,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)(Cds-Cd)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Cd 0 {1,S}
        5     Cb 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Cb",
    index=737,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cb 0 {1,S}
        6     Cd 0 {3,D}
        7     Cd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=738,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)(Cds-Cds)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cb 0 {1,S}
        6     Cdd 0 {3,D}
        7     Cd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cds)Cb",
    index=739,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cb 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cd 0 {4,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)",
    index=740,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cds)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cb 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cd 0 {4,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Cb",
    index=741,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)(Cds-Cdd)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cb 0 {1,S}
        6     Cdd 0 {3,D}
        7     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb",
    index=742,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cb 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    index=743,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cb 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)Cb",
    index=744,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Cb 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Cb",
    index=745,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)(Cds-Cd)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Cd 0 {1,S}
        5     Cb 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cb",
    index=746,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cb 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cd 0 {4,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=747,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cb 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)Cb",
    index=748,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cb 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cdd 0 {4,D} {9,D}
        9     Od 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)",
    index=749,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cb 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cdd 0 {4,D} {9,D}
        9     C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cb",
    index=750,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd)(Cds-Cdd)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cb 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D}
        8     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb",
    index=751,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cb 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D} {9,D}
        8     Cdd 0 {4,D} {10,D}
        9     Od 0 {7,D}
        10    Od 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    index=752,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Cd)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cb 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D} {9,D}
        8     Cdd 0 {4,D} {10,D}
        9     Od 0 {7,D}
        10    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Cb",
    index=753,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cb 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D} {9,D}
        8     Cdd 0 {4,D} {10,D}
        9     C 0 {7,D}
        10    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cb",
    index=754,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cb 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        8     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb",
    index=755,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cb 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     Od 0 {6,D}
        10    Od 0 {7,D}
        11    Od 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    index=756,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cb 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     Od 0 {6,D}
        10    Od 0 {7,D}
        11    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)Cb",
    index=757,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cb 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     Od 0 {6,D}
        10    C 0 {7,D}
        11    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Cb",
    index=758,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Cb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Cb 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     C 0 {6,D}
        10    C 0 {7,D}
        11    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Cb",
    index=759,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCtCdsCds",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     {Cd,CO} 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtCt",
    index=760,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)CtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)(Cds-Cds)",
    index=761,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)CtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CtCt",
    index=762,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)CtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=763,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)CtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)CtCt",
    index=764,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)CtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)",
    index=765,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)CtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CtCt",
    index=766,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)CtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtCt",
    index=767,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)CtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.61,7.3,8.97,9.69,9.84,9.42,7.36],"cal/(mol*K)"),
        H298=(5.48,"kcal/mol"),
        S298=(-34.5,"cal/(mol*K)"),
    ),
    index=768,
    short_comment="Cs-CtCtCdCd BOZZELLI =3D Cs/Cs/Cd/Ct2 + (Cs/Cs3/Cd - Cs/Cs4)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cds)CtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cds)CtCt",
    index=769,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cds)CtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)",
    index=770,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cds)CtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtCt",
    index=771,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)CtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtCt",
    index=772,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    index=773,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CtCt",
    index=774,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtCt",
    index=775,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCtCdsCds",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     {Cd,CO} 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbCt",
    index=776,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)CbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Ct",
    index=777,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)CbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CbCt",
    index=778,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)CbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        6     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Ct",
    index=779,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)CbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)CbCt",
    index=780,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)CbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)Ct",
    index=781,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)CbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CbCt",
    index=782,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)CbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbCt",
    index=783,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)CbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.61,7.3,8.97,9.69,9.84,9.42,7.36],"cal/(mol*K)"),
        H298=(5.48,"kcal/mol"),
        S298=(-34.5,"cal/(mol*K)"),
    ),
    index=784,
    short_comment="Cs-CbCtCdCd BOZZELLI =3D Cs/Cs/Cb/Cd2 + (Cs/Cs3/Ct - Cs/Cs4)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cds)CbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cds)CbCt",
    index=785,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cds)CbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Ct",
    index=786,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cds)CbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbCt",
    index=787,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)CbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbCt",
    index=788,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)Ct",
    index=789,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CbCt",
    index=790,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbCt",
    index=791,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbCdsCds",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     {Cd,CO} 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbCb",
    index=792,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)CbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)(Cds-Cds)",
    index=793,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)CbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CbCb",
    index=794,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)CbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        6     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=795,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)CbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        6     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)CbCb",
    index=796,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)CbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)",
    index=797,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)CbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CbCb",
    index=798,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)CbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbCb",
    index=799,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)CbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.61,7.3,8.97,9.69,9.84,9.42,7.36],"cal/(mol*K)"),
        H298=(5.48,"kcal/mol"),
        S298=(-34.5,"cal/(mol*K)"),
    ),
    index=800,
    short_comment="Cs-CbCbCdCd BOZZELLI =3D Cs/Cs/Cb2/Cd + (Cs/Cs3/Cd - Cs/Cs4)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cds)CbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cds)CbCb",
    index=801,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cds)CbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)",
    index=802,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cds)CbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbCb",
    index=803,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)CbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbCb",
    index=804,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)",
    index=805,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CbCb",
    index=806,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbCb",
    index=807,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCtCtCds",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     {Cd,CO} 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtCtCt",
    index=808,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CtCtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=809,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CtCtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-CtCtCtCds",
    index=0,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CtCtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cd 0 {2,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=810,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CtCtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CtCtCt",
    index=811,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CtCtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)",
    index=812,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CtCtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CtCtCt",
    index=813,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCtCtCds",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     {Cd,CO} 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCtCt",
    index=814,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CbCtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CtCt",
    index=815,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CbCtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCtCt",
    index=816,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CbCtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cd 0 {2,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtCt",
    index=817,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CbCtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CbCtCt",
    index=818,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CbCtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CtCt",
    index=819,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CbCtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CbCtCt",
    index=820,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbCtCds",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     {Cd,CO} 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCbCt",
    index=821,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CbCbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Ct",
    index=822,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CbCbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCbCt",
    index=823,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CbCbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        6     Cd 0 {2,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Ct",
    index=824,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CbCbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CbCbCt",
    index=825,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CbCbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Ct",
    index=826,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CbCbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CbCbCt",
    index=827,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbCbCds",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     {Cd,CO} 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCbCb",
    index=828,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CbCbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=829,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CbCbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCbCb",
    index=830,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CbCbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        6     Cd 0 {2,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=831,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CbCbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CbCbCb",
    index=832,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CbCbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)",
    index=833,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CbCbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CbCbCb",
    index=834,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCtCtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=835,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCtCtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtCtCt",
    index=836,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbCtCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtCt",
    index=837,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbCbCt",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Ct 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Ct",
    index=838,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbCbCb",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Cb 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)(Cds-Cds)",
    index=839,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CCCOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     C 0 {1,S}
        3     C 0 {1,S}
        4     C 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-CsCsCsOs",
    index=840,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CsCsCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cs 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.33,6.19,7.25,7.7,8.2,8.24,8.24],"cal/(mol*K)"),
        H298=(-6.6,"kcal/mol"),
        S298=(-33.56,"cal/(mol*K)"),
    ),
    index=841,
    short_comment="Cs-OCsCsCs BENSON !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsCsCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     {Cd,CO} 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CsCsOs",
    index=842,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CsCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.99,6.04,7.43,8.26,8.92,8.96,8.23],"cal/(mol*K)"),
        H298=(-3.6,"kcal/mol"),
        S298=(-34.72,"cal/(mol*K)"),
    ),
    index=843,
    short_comment="Cs-OCOCsCs Hf BENSON S,Cp =3D C/Cd/C3",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CsCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CsCsOs",
    index=844,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CsCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.63,6.79,7.95,8.4,8.8,8.44,8.44],"cal/(mol*K)"),
        H298=(-6.6,"kcal/mol"),
        S298=(-32.56,"cal/(mol*K)"),
    ),
    index=845,
    short_comment="Cs-OCdCsCs BOZZELLI C/C3/O - (C/C3/H - C/Cb/C2/H), Hf-1 !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CsCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CsCsOs",
    index=846,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CsCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([8.39,9.66,10.03,10.07,9.64,9.26,8.74],"cal/(mol*K)"),
        H298=(-9.725,"kcal/mol"),
        S298=(-36.5,"cal/(mol*K)"),
    ),
    index=847,
    short_comment="{C/CCO/O/C2} RAMAN & GREEN JPCA 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CsCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CsCsOs",
    index=848,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-OsCtCsCs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Os 0 {1,S}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     Cs 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CsCsOs",
    index=849,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCsCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Cs 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.63,6.79,7.95,8.4,8.8,8.44,8.44],"cal/(mol*K)"),
        H298=(-6.6,"kcal/mol"),
        S298=(-32.56,"cal/(mol*K)"),
    ),
    index=850,
    short_comment="Cs-OCbCsCs BOZZELLI C/C3/O - (C/C3/H - C/Cb/C2/H), Hf-1 !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsCdsCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     {Cd,CO} 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CsOs",
    index=851,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)CsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-CsCsCsOs",
    index=852,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)CsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CsOs",
    index=853,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)CsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Od)CsCsOs",
    index=854,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)CsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)CsOs",
    index=855,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)CsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsCsOs",
    index=856,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)CsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CsOs",
    index=857,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)CsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CsOs",
    index=858,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)CsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.61,5.98,7.51,8.37,9,9.02,8.34],"cal/(mol*K)"),
        H298=(-8.01,"kcal/mol"),
        S298=(-34.34,"cal/(mol*K)"),
    ),
    index=859,
    short_comment="Cs-OCdCdCs Hf jwb 697 S,Cp from C/Cd2/C2",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cds)CsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cds)CsOs",
    index=860,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cds)CsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsCsOs",
    index=861,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cds)CsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CsOs",
    index=862,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)CsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CsOs",
    index=863,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CsOs",
    index=864,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CsOs",
    index=865,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CsOs",
    index=866,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCdsCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtCsOs",
    index=867,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CtCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CsOs",
    index=868,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CtCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtCsOs",
    index=869,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CtCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CsOs",
    index=870,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CtCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CtCsOs",
    index=871,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CtCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CsOs",
    index=872,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CtCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CtCsOs",
    index=873,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCdsCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCsOs",
    index=874,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CbCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CsOs",
    index=875,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CbCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCsOs",
    index=876,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CbCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CsOs",
    index=877,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CbCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CbCsOs",
    index=878,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CbCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CsOs",
    index=879,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CbCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CbCsOs",
    index=880,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCtCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CsOs",
    index=881,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCtCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtCsOs",
    index=882,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbCsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Cs 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CsOs",
    index=883,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsCdsCdsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     {Cd,CO} 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Os",
    index=884,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Od)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     CO 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-CsCsCsOs",
    index=885,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cd)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Os",
    index=886,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Os 0 {1,S}
        6     Cd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)CsOs",
    index=887,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Os 0 {1,S}
        6     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)Os",
    index=888,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Od)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Os 0 {1,S}
        6     Cdd 0 {4,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsCsOs",
    index=889,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)(Cds-Cdd-Cd)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cd 0 {1,S} {6,D}
        5     Os 0 {1,S}
        6     Cdd 0 {4,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Os",
    index=890,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)(Cds-Cd)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Cd 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Os",
    index=891,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Os 0 {1,S}
        6     Cd 0 {3,D}
        7     Cd 0 {4,D}
        """,
    node="Cs-(Cds-Od)CsCsOs",
    index=892,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)(Cds-Cds)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Os 0 {1,S}
        6     Cdd 0 {3,D}
        7     Cd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cds)Os",
    index=893,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Os 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cd 0 {4,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Od)CsOs",
    index=894,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cds)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Os 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cd 0 {4,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Os",
    index=895,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)(Cds-Cdd)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Os 0 {1,S}
        6     Cdd 0 {3,D}
        7     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Os",
    index=896,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Os 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CsOs",
    index=897,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Os 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)Os",
    index=898,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cd 0 {1,S} {7,D}
        5     Os 0 {1,S}
        6     Cdd 0 {3,D} {8,D}
        7     Cdd 0 {4,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Os",
    index=899,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)(Cds-Cd)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Cd 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Os",
    index=900,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cd 0 {4,D}
        """,
    node="Cs-CsCsCsOs",
    index=901,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)Os",
    index=902,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cdd 0 {4,D} {9,D}
        9     Od 0 {8,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsCsOs",
    index=903,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Cd)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        8     Cdd 0 {4,D} {9,D}
        9     C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Os",
    index=904,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd)(Cds-Cdd)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D}
        8     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Os",
    index=905,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D} {9,D}
        8     Cdd 0 {4,D} {10,D}
        9     Od 0 {7,D}
        10    Od 0 {8,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CsOs",
    index=906,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Cd)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D} {9,D}
        8     Cdd 0 {4,D} {10,D}
        9     Od 0 {7,D}
        10    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Os",
    index=907,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        7     Cdd 0 {3,D} {9,D}
        8     Cdd 0 {4,D} {10,D}
        9     C 0 {7,D}
        10    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Os",
    index=908,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)(Cds-Cdd)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        8     Cdd 0 {4,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Os",
    index=909,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Od)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     Od 0 {6,D}
        10    Od 0 {7,D}
        11    Od 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Os",
    index=910,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)(Cds-Cdd-Cd)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     Od 0 {6,D}
        10    Od 0 {7,D}
        11    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)Os",
    index=911,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     Od 0 {6,D}
        10    C 0 {7,D}
        11    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Os",
    index=912,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)(Cds-Cdd-Cd)Os",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cd 0 {1,S} {8,D}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {9,D}
        7     Cdd 0 {3,D} {10,D}
        8     Cdd 0 {4,D} {11,D}
        9     C 0 {6,D}
        10    C 0 {7,D}
        11    C 0 {8,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Os",
    index=913,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCdsCdsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtOs",
    index=914,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)CtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Os",
    index=915,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)CtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CtOs",
    index=916,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)CtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Os",
    index=917,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)CtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)CtOs",
    index=918,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)CtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)Os",
    index=919,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)CtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CtOs",
    index=920,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)CtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtOs",
    index=921,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)CtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Os",
    index=922,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cds)CtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cds)CtOs",
    index=923,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cds)CtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Os",
    index=924,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cds)CtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtOs",
    index=925,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)CtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtOs",
    index=926,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)Os",
    index=927,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CtOs",
    index=928,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtOs",
    index=929,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCdsCdsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbOs",
    index=930,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)CbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Od)(Cds-Cds)Os",
    index=931,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)CbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CbOs",
    index=932,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)CbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        6     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Os",
    index=933,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)CbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)CbOs",
    index=934,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)CbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Od)(Cds-Cds)Os",
    index=935,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)CbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CbOs",
    index=936,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)CbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbOs",
    index=937,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)CbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Os",
    index=938,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cds)CbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cds)CbOs",
    index=939,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cds)CbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Os",
    index=940,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cds)CbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbOs",
    index=941,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)CbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbOs",
    index=942,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)CbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cdd-Od)(Cds-Cdd-Od)Os",
    index=943,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)CbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CbOs",
    index=944,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)CbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CbOs",
    index=945,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCtCdsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtCtOs",
    index=946,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CtCtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Os",
    index=947,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CtCtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtCtOs",
    index=948,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CtCtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Os",
    index=949,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CtCtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CtCtOs",
    index=950,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CtCtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Os",
    index=951,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CtCtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CtCtOs",
    index=952,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCtCdsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCtOs",
    index=953,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CbCtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)CtOs",
    index=954,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CbCtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCtOs",
    index=955,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CbCtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtOs",
    index=956,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CbCtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CbCtOs",
    index=957,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CbCtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)CtOs",
    index=958,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CbCtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CbCtOs",
    index=959,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbCdsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     {Cd,CO} 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCbOs",
    index=960,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CbCbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)(Cds-Cds)Os",
    index=961,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CbCbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbCbOs",
    index=962,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CbCbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Os",
    index=963,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CbCbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CbCbOs",
    index=964,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CbCbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cdd-Od)Os",
    index=965,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CbCbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CbCbOs",
    index=966,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCtCtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Os",
    index=967,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCtCtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtCtOs",
    index=968,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbCtOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Ct 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)CtOs",
    index=969,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbCbOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Cb 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)(Cds-Cds)Os",
    index=970,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CCOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     C 0 {1,S}
        3     C 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-CsCsOsOs",
    index=971,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CsCsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cs 0 {1,S}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.8,6.09,7.3,7.78,8.24,8.24,8.24],"cal/(mol*K)"),
        H298=(-9.77,"kcal/mol"),
        S298=(-33.18,"cal/(mol*K)"),
    ),
    index=972,
    short_comment="Cs-OOCsCs BOZZELLI =3D C/C3/O - (C/C2/H2 - C/C/H2/O) !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsCsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     {Cd,CO} 0 {1,S}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CsOsOs",
    index=973,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-CsCsOsOs",
    index=974,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CsOsOs",
    index=975,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        """,
    node="Cs-CsCsOsOs",
    index=976,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CsOsOs",
    index=977,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CsOsOs",
    index=978,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CsOsOs",
    index=979,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsCdsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     {Cd,CO} 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)OsOs",
    index=980,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)OsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-CsCsOsOs",
    index=981,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)OsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)OsOs",
    index=982,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)OsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Od)CsOsOs",
    index=983,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)OsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)OsOs",
    index=984,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)OsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsOsOs",
    index=985,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)OsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)OsOs",
    index=986,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)OsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)OsOs",
    index=987,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)OsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-CsCsOsOs",
    index=988,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cds)OsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cds)OsOs",
    index=989,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cds)OsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsOsOs",
    index=990,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cds)OsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)OsOs",
    index=991,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)OsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)OsOs",
    index=992,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)OsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)OsOs",
    index=993,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)OsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)OsOs",
    index=994,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)OsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)OsOs",
    index=995,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Ct 0 {1,S}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CsOsOs",
    index=996,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCdsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Ct 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtOsOs",
    index=997,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CtOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)OsOs",
    index=998,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CtOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtOsOs",
    index=999,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CtOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)OsOs",
    index=1000,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CtOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CtOsOs",
    index=1001,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CtOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)OsOs",
    index=1002,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CtOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CtOsOs",
    index=1003,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCtOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)OsOs",
    index=1004,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CsOsOs",
    index=1005,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCdsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbOsOs",
    index=1006,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CbOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cb 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)OsOs",
    index=1007,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CbOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cb 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbOsOs",
    index=1008,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CbOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)OsOs",
    index=1009,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CbOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CbOsOs",
    index=1010,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CbOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)OsOs",
    index=1011,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CbOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CbOsOs",
    index=1012,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCtOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtOsOs",
    index=1013,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)OsOs",
    index=1014,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-COsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     C 0 {1,S}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-CsOsOsOs",
    index=1015,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CsOsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cs 0 {1,S}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.33,6.19,7.25,7.7,8.2,8.24,8.24],"cal/(mol*K)"),
        H298=(-19,"kcal/mol"),
        S298=(-33.56,"cal/(mol*K)"),
    ),
    index=1016,
    short_comment="Cs-OOOCs BOZZELLI est !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsOsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     {Cd,CO} 0 {1,S}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)OsOsOs",
    index=1017,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)OsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-CsOsOsOs",
    index=1018,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)OsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)OsOsOs",
    index=1019,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)OsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cd 0 {2,D}
        """,
    node="Cs-CsOsOsOs",
    index=1020,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)OsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)OsOsOs",
    index=1021,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)OsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)OsOsOs",
    index=1022,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)OsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)OsOsOs",
    index=1023,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtOsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Ct 0 {1,S}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)OsOsOs",
    index=1024,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbOsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    node="Cs-(Cds-Cds)OsOsOs",
    index=1025,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-OsOsOsOs",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Os 0 {1,S}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     Os 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.33,6.13,7.25,7.7,8.2,8.24,8.24],"cal/(mol*K)"),
        H298=(-23,"kcal/mol"),
        S298=(-35.56,"cal/(mol*K)"),
    ),
    index=1026,
    short_comment="Cs-OOOO BOZZELLI est !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-COsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     C 0 {1,S}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-CsOsOsH",
    index=1027,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CsOsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cs 0 {1,S}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.25,7.1,8.81,9.55,10.31,11.05,11.05],"cal/(mol*K)"),
        H298=(-16,"kcal/mol"),
        S298=(-12.07,"cal/(mol*K)"),
    ),
    index=1028,
    short_comment="Cs-OOCsH BENSON Hf, BOZZELLI C/C3/H - C/C2/O/H !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsOsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     {Cd,CO} 0 {1,S}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)OsOsH",
    index=1029,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)OsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-CsOsOsH",
    index=1030,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)OsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)OsOsH",
    index=1031,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)OsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        """,
    node="Cs-CsOsOsH",
    index=1032,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)OsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)OsOsH",
    index=1033,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)OsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cds)OsOsH",
    index=1034,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)OsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)OsOsH",
    index=1035,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtOsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)OsOsH",
    index=1036,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbOsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Os 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)OsOsH",
    index=1037,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CCOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     C 0 {1,S}
        3     C 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-CsCsOsH",
    index=1038,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CsCsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cs 0 {1,S}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.8,6.64,8.1,8.73,9.81,10.4,11.51],"cal/(mol*K)"),
        H298=(-7.2,"kcal/mol"),
        S298=(-11,"cal/(mol*K)"),
    ),
    index=1039,
    short_comment="Cs-OCsCs BENSON: Cp1500 =3D Cp1000*(Cp1500/Cp1000: C/C2Cd/H)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsCsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     {Cd,CO} 0 {1,S}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CsOsH",
    index=1040,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.47,6.82,8.45,9.17,10.24,10.8,11.02],"cal/(mol*K)"),
        H298=(-6,"kcal/mol"),
        S298=(-11.1,"cal/(mol*K)"),
    ),
    index=1041,
    short_comment="Cs-OCOCsH BOZZELLI",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CsOsH",
    index=1042,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.47,6.82,8.45,9.17,10.24,10.8,11.02],"cal/(mol*K)"),
        H298=(-6,"kcal/mol"),
        S298=(-11.1,"cal/(mol*K)"),
    ),
    index=1043,
    short_comment="Cs-OCdCsH BOZZELLI",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CsOsH",
    index=1044,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([7.2,8.49,9.33,9.92,10.5,10.92,11.71],"cal/(mol*K)"),
        H298=(-8.37,"kcal/mol"),
        S298=(-13.04,"cal/(mol*K)"),
    ),
    index=1045,
    short_comment="{C/CCO/O/C/H} RAMAN & GREEN JPCA 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CsOsH",
    index=1046,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsCdsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     {Cd,CO} 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)OsH",
    index=1047,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Od)OsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-CsCsOsH",
    index=1048,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cd)OsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)OsH",
    index=1049,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cds)OsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Od)CsOsH",
    index=1050,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd)OsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cdd-Cd)OsH",
    index=1051,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Od)OsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsOsH",
    index=1052,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)(Cds-Cdd-Cd)OsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S} {6,D}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {3,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)OsH",
    index=1053,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)(Cds-Cd)OsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)OsH",
    index=1054,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)(Cds-Cds)OsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.21,6.6,8.26,9.05,10.23,10.86,11.04],"cal/(mol*K)"),
        H298=(-6.67,"kcal/mol"),
        S298=(-10.42,"cal/(mol*K)"),
    ),
    index=1055,
    short_comment="Cs-OCdCdH BOZZELLI",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cds)OsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cds)OsH",
    index=1056,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cds)OsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)CsOsH",
    index=1057,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cds)OsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cd 0 {3,D}
        8     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)OsH",
    index=1058,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)(Cds-Cdd)OsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        7     Cdd 0 {3,D}
        """,
    node="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)OsH",
    index=1059,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Od)OsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     Od 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)OsH",
    index=1060,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)(Cds-Cdd-Cd)OsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     Od 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)OsH",
    index=1061,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)(Cds-Cdd-Cd)OsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cd 0 {1,S} {7,D}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {8,D}
        7     Cdd 0 {3,D} {9,D}
        8     C 0 {6,D}
        9     C 0 {7,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)OsH",
    index=1062,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CsOsH",
    index=1063,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCdsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtOsH",
    index=1064,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CtOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)OsH",
    index=1065,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CtOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtOsH",
    index=1066,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CtOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)OsH",
    index=1067,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CtOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CtOsH",
    index=1068,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CtOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)OsH",
    index=1069,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CtOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CtOsH",
    index=1070,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtCtOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)OsH",
    index=1071,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Cs 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.47,6.82,8.45,9.17,10.24,10.8,11.02],"cal/(mol*K)"),
        H298=(-6,"kcal/mol"),
        S298=(-11.1,"cal/(mol*K)"),
    ),
    index=1072,
    short_comment="Cs-OCbCsH BOZZELLI =3D C/Cd/C/H/O Jul 91",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCdsOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbOsH",
    index=1073,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)CbOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Cb 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Od)(Cds-Cds)OsH",
    index=1074,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)CbOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S}
        3     Cb 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CbOsH",
    index=1075,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)CbOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)OsH",
    index=1076,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)CbOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)CbOsH",
    index=1077,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)CbOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    node="Cs-(Cds-Cdd-Od)(Cds-Cds)OsH",
    index=1078,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)CbOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Cb 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)CbOsH",
    index=1079,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCtOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Ct 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)CtOsH",
    index=1080,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbCbOsH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        4     Os 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)(Cds-Cds)OsH",
    index=1081,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-COsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     C 0 {1,S}
        3     Os 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-CsOsHH",
    index=1082,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CsOsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cs 0 {1,S}
        3     Os 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.99,6.85,8.3,9.43,11.11,12.33,12.33],"cal/(mol*K)"),
        H298=(-8.1,"kcal/mol"),
        S298=(9.8,"cal/(mol*K)"),
    ),
    index=1083,
    short_comment="Cs-OCsHH BENSON !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CdsOsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     {Cd,CO} 0 {1,S}
        3     Os 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)OsHH",
    index=1084,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Od)OsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     CO 0 {1,S}
        3     Os 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([2.95,6.74,8.53,9.8,11.61,12.65,14.4],"cal/(mol*K)"),
        H298=(-5.28,"kcal/mol"),
        S298=(10.17,"cal/(mol*K)"),
    ),
    index=1085,
    short_comment="Cs-OCOHH jwl 99 cdsQ cc*ocq",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cd)OsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S}
        3     Os 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)OsHH",
    index=1086,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cds)OsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Os 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cd 0 {2,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.12,6.86,8.32,9.49,11.22,12.48,14.4],"cal/(mol*K)"),
        H298=(-6.76,"kcal/mol"),
        S298=(9.8,"cal/(mol*K)"),
    ),
    index=1087,
    short_comment="Cs-OCdHH BOZZELLI Hf PEDLEY c*ccoh C/C/Cd/H2",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd)OsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Os 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D}
        """,
    node="Cs-(Cds-Cdd-Cd)OsHH",
    index=1088,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Od)OsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cd 0 {1,S} {6,D}
        3     Os 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     Od 0 {6,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([7.15,8.67,9.75,10.65,11.93,12.97,14.86],"cal/(mol*K)"),
        H298=(-8.68,"kcal/mol"),
        S298=(8.43,"cal/(mol*K)"),
    ),
    index=1089,
    short_comment="{C/CCO/O/H2} RAMAN & GREEN JPCA 2002, 106, 7937-7949",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-(Cds-Cdd-Cd)OsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {5,S} {4,S}
        2     Cd 0 {1,S} {6,D}
        3     Os 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        6     Cdd 0 {2,D} {7,D}
        7     C 0 {6,D}
        """,
    node="Cs-(Cds-Cds)OsHH",
    index=1090,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CtOsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Ct 0 {1,S}
        3     Os 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.12,6.86,8.32,9.49,11.22,12.48,14.4],"cal/(mol*K)"),
        H298=(-6.76,"kcal/mol"),
        S298=(9.8,"cal/(mol*K)"),
    ),
    index=1091,
    short_comment="Cs-OCtHH BOZZELLI assigned C/Cd/H2/O",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Cs-CbOsHH",
    group=
        """
        1  *  Cs 0 {2,S} {3,S} {4,S} {5,S}
        2     Cb 0 {1,S}
        3     Os 0 {1,S}
        4     H 0 {1,S}
        5     H 0 {1,S}
        """,
    node="Cs-(Cds-Cds)OsHH",
    index=1092,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="O",
    group=
        """
        1  *  O 0
        """,
    node="Os-CsCs",
    index=1093,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Od",
    group=
        """
        1  *  Od 0
        """,
    node="Od-Cd",
    index=1094,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Od-Cd",
    group=
        """
        1  *  Od 0 {2,D}
        2     {Cd,CO} 0 {1,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([0,0,0,0,0,0,0],"cal/(mol*K)"),
        H298=(0,"kcal/mol"),
        S298=(0,"cal/(mol*K)"),
    ),
    index=1095,
    short_comment="In this case the C is treated as the central atom",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Od-Od",
    group=
        """
        1  *  Od 0 {2,D}
        2     Od 0 {1,D}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.5,3.61,3.72,3.825,4.035,4.175,4.36],"cal/(mol*K)"),
        H298=(0,"kcal/mol"),
        S298=(25.19,"cal/(mol*K)"),
    ),
    index=1096,
    short_comment="O2 Kee et al., SAND87-8215B, 1994, we cut it half to account for adding two single Od groups, also the symmetric number is considered too. S(group) = (S(o2) + Rln(sigma))/2",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os",
    group=
        """
        1  *  Os 0
        """,
    node="Os-(Cds-Cd)(Cds-Cd)",
    index=1097,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-HH",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     H 0 {1,S}
        3     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([8.03,8.19,8.42,8.68,9.26,9.86,11.26],"cal/(mol*K)"),
        H298=(-57.8,"kcal/mol"),
        S298=(46.51,"cal/(mol*K)"),
    ),
    index=1098,
    short_comment="O-HH WATER. !!!Using NIST value for H2O, S(group) = S(H2O) + Rln(2)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-OsH",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Os 0 {1,S}
        3     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([5.21,5.72,6.17,6.66,7.15,7.61,8.43],"cal/(mol*K)"),
        H298=(-16.3,"kcal/mol"),
        S298=(27.83,"cal/(mol*K)"),
    ),
    index=1099,
    short_comment="O-OH SANDIA 1/2*H2O2",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-OsOs",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Os 0 {1,S}
        3     Os 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([2.2,3.64,4.2,4.34,4.62,4.9,4.9],"cal/(mol*K)"),
        H298=(8.85,"kcal/mol"),
        S298=(9.4,"cal/(mol*K)"),
    ),
    index=1100,
    short_comment="O-OO LAY 1997=20 !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-CH",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     C 0 {1,S}
        3     H 0 {1,S}
        """,
    node="Os-CsH",
    index=1101,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-CtH",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Ct 0 {1,S}
        3     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.3,4.5,4.82,5.23,6.02,6.61,7.44],"cal/(mol*K)"),
        H298=(-37.9,"kcal/mol"),
        S298=(29.1,"cal/(mol*K)"),
    ),
    index=1102,
    short_comment="O-CtH BENSON (Assigned O-CsH)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-CdsH",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     {Cd,CO} 0 {1,S}
        3     H 0 {1,S}
        """,
    node="Os-(Cds-Cd)H",
    index=1103,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-(Cds-Od)H",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     CO 0 {1,S}
        3     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.8,5,5.8,6.3,7.2,7.8,7.8],"cal/(mol*K)"),
        H298=(-58.1,"kcal/mol"),
        S298=(24.5,"cal/(mol*K)"),
    ),
    index=1104,
    short_comment="O-COH !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-(Cds-Cd)H",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Cd 0 {1,S}
        3     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.3,4.5,4.82,5.23,6.02,6.61,7.44],"cal/(mol*K)"),
        H298=(-37.9,"kcal/mol"),
        S298=(29.1,"cal/(mol*K)"),
    ),
    index=1105,
    short_comment="O-CdH BENSON (Assigned O-CsH)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-CsH",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Cs 0 {1,S}
        3     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.3,4.5,4.82,5.23,6.02,6.61,7.44],"cal/(mol*K)"),
        H298=(-37.9,"kcal/mol"),
        S298=(29.07,"cal/(mol*K)"),
    ),
    index=1106,
    short_comment="O-CsH BENSON",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-CbH",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Cb 0 {1,S}
        3     H 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([4.3,4.5,4.82,5.23,6.02,6.61,7.44],"cal/(mol*K)"),
        H298=(-37.9,"kcal/mol"),
        S298=(29.1,"cal/(mol*K)"),
    ),
    index=1107,
    short_comment="O-CbH BENSON (Assigned O-CsH)",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-OsC",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Os 0 {1,S}
        3     C 0 {1,S}
        """,
    node="Os-OsCs",
    index=1108,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-OsCt",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Os 0 {1,S}
        3     Ct 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.9,4.31,4.6,4.84,5.32,5.8,5.8],"cal/(mol*K)"),
        H298=(7,"kcal/mol"),
        S298=(10.8,"cal/(mol*K)"),
    ),
    index=1109,
    short_comment="O-OCb Hf JWB plot S,Cp assigned O/O/Cd !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-OsCds",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Os 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        """,
    node="Os-Os(Cds-Cd)",
    index=1110,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-Os(Cds-Od)",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Os 0 {1,S}
        3     CO 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.53,5.02,5.79,6.08,6.54,6.49,6.49],"cal/(mol*K)"),
        H298=(-23.22,"kcal/mol"),
        S298=(9.11,"cal/(mol*K)"),
    ),
    index=1111,
    short_comment="O-OCO jwl cbsQ 99 cqcho=20 !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-Os(Cds-Cd)",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Os 0 {1,S}
        3     Cd 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.5,3.87,3.95,4.15,4.73,4.89,4.89],"cal/(mol*K)"),
        H298=(1.64,"kcal/mol"),
        S298=(10.12,"cal/(mol*K)"),
    ),
    index=1112,
    short_comment="O-OCd WESTMORELAND S,Cp LAY'9405 !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-OsCs",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Os 0 {1,S}
        3     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.9,4.31,4.6,4.84,5.32,5.8,5.8],"cal/(mol*K)"),
        H298=(-5.4,"kcal/mol"),
        S298=(8.54,"cal/(mol*K)"),
    ),
    index=1113,
    short_comment="O-OCs LAY 1997 !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-OsCb",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Os 0 {1,S}
        3     Cb 0 {1,S}
        """,
    node="Os-Os(Cds-Cd)",
    index=1114,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-CC",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     C 0 {1,S}
        3     C 0 {1,S}
        """,
    node="Os-(Cds-Cd)(Cds-Cd)",
    index=1115,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-CtCt",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Ct 0 {1,S}
        3     Ct 0 {1,S}
        """,
    node="Os-(Cds-Cd)(Cds-Cd)",
    index=1116,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-CtCds",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Ct 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        """,
    node="Os-(Cds-Cd)(Cds-Cd)",
    index=1117,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-Ct(Cds-Od)",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Ct 0 {1,S}
        3     CO 0 {1,S}
        """,
    node="Os-(Cds-Cd)(Cds-Cd)",
    index=1118,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-Ct(Cds-Cd)",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Ct 0 {1,S}
        3     Cd 0 {1,S}
        """,
    node="Os-(Cds-Cd)(Cds-Cd)",
    index=1119,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-CtCs",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Ct 0 {1,S}
        3     Cs 0 {1,S}
        """,
    node="Os-Cs(Cds-Cd)",
    index=1120,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-CtCb",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Ct 0 {1,S}
        3     Cb 0 {1,S}
        """,
    node="Os-(Cds-Cd)(Cds-Cd)",
    index=1121,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-CdsCds",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     {Cd,CO} 0 {1,S}
        3     {Cd,CO} 0 {1,S}
        """,
    node="Os-(Cds-Cd)(Cds-Cd)",
    index=1122,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-(Cds-Od)(Cds-Od)",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     CO 0 {1,S}
        3     CO 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.45,4.74,5.28,5.74,5.89,6.1,6.1],"cal/(mol*K)"),
        H298=(-46,"kcal/mol"),
        S298=(2.5,"cal/(mol*K)"),
    ),
    index=1123,
    short_comment="O-COCO Hf BENSON S,Cp Mopac=3D",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-(Cds-Od)(Cds-Cd)",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     CO 0 {1,S}
        3     Cd 0 {1,S}
        """,
    node="Os-(Cds-Cd)(Cds-Cd)",
    index=1124,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-(Cds-Cd)(Cds-Cd)",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Cd 0 {1,S}
        3     Cd 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.4,3.7,3.7,3.8,4.4,4.6,4.8],"cal/(mol*K)"),
        H298=(-19.61,"kcal/mol"),
        S298=(10,"cal/(mol*K)"),
    ),
    index=1125,
    short_comment="O-CdCd BOZZELLI",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-CdsCs",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     {Cd,CO} 0 {1,S}
        3     Cs 0 {1,S}
        """,
    node="Os-Cs(Cds-Cd)",
    index=1126,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-Cs(Cds-Od)",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Cs 0 {1,S}
        3     CO 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.91,4.31,4.6,4.84,5.32,5.8,5.8],"cal/(mol*K)"),
        H298=(-42.19,"kcal/mol"),
        S298=(8.4,"cal/(mol*K)"),
    ),
    index=1127,
    short_comment="O-COCs BOZZELLI Jul91 S,Cp ABaldwin O/Cs/O !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-Cs(Cds-Cd)",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Cs 0 {1,S}
        3     Cd 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.91,4.31,4.6,4.84,5.32,5.8,5.8],"cal/(mol*K)"),
        H298=(-23.73,"kcal/mol"),
        S298=(9.7,"cal/(mol*K)"),
    ),
    index=1128,
    short_comment="O-CdCs Hf RADOM vin-oh S A.Baldwin O/Cs/O !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-CdsCb",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     {Cd,CO} 0 {1,S}
        3     Cb 0 {1,S}
        """,
    node="Os-CC",
    index=0,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-Cb(Cds-Od)",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Cb 0 {1,S}
        3     CO 0 {1,S}
        """,
    node="Os-CdsCb",
    index=0,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-Cb(Cds-Cd)",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Cb 0 {1,S}
        3     Cd 0 {1,S}
        """,
    node="Os-CdsCb",
    index=0,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-CsCs",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Cs 0 {1,S}
        3     Cs 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.4,3.7,3.7,3.8,4.4,4.6,4.6],"cal/(mol*K)"),
        H298=(-23.2,"kcal/mol"),
        S298=(8.68,"cal/(mol*K)"),
    ),
    index=1129,
    short_comment="O-CsCs BENSON !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-CsCb",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Cs 0 {1,S}
        3     Cb 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([3.4,3.7,3.7,3.8,4.4,4.6,4.6],"cal/(mol*K)"),
        H298=(-22.6,"kcal/mol"),
        S298=(9.7,"cal/(mol*K)"),
    ),
    index=1130,
    short_comment="O-CbCs REID, PRAUSNITZ and SHERWOOD !!!WARNING! Cp1500 value taken as Cp1000",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Os-CbCb",
    group=
        """
        1  *  Os 0 {2,S} {3,S}
        2     Cb 0 {1,S}
        3     Cb 0 {1,S}
        """,
    model=ThermoGAModel(
        Tdata=([300,400,500,600,800,1000,1500],"K"),
        Cpdata=([1.19,-0.24,-0.72,-0.51,0.43,1.36,1.75],"cal/(mol*K)"),
        H298=(-18.77,"kcal/mol"),
        S298=(13.59,"cal/(mol*K)"),
    ),
    index=1131,
    short_comment="O-CbCb CHERN 1/97 Hf PEDLEY, Mopac",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="Si",
    group=
        """
        1  *  Si 0
        """,
    node="Cs-HHHH",
    index=1132,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

thermo(
    label="S",
    group=
        """
        1  *  S 0
        """,
    node="Os-HH",
    index=1133,
    short_comment="",
    long_comment=
        """
        """,
    history=[
    ],
)

