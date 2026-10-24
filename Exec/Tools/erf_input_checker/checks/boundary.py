"""Boundary conditions and the surface layer."""

from ..checker import Checker, register

FACE_TYPES = {"inflow", "outflow", "inflow_outflow", "open", "ho_outflow", "slipwall",
              "noslipwall", "symmetry", "surface_layer", "most"}
FACES = ("xlo", "xhi", "ylo", "yhi", "zlo", "zhi")
COLUMN_PBL = {"mynn25", "mynnedmf", "myj", "ysu", "ysunew", "mrf", "native_shoc", "eamxx_shoc", "shoc"}


@register
class BoundaryChecker(Checker):
    name = "boundary"

    def fkey(self, face: str, suffix: str):
        """The key as the deck spells it: ERF reads erf.<face>.<suffix> first, then <face>.<suffix>."""
        for k in ("erf.%s.%s" % (face, suffix), "%s.%s" % (face, suffix)):
            if self.deck.has(k):
                return k
        return None

    def fstr(self, face: str, suffix: str, default=None):
        k = self.fkey(face, suffix)
        return self.deck.get_str(k, default) if k else default

    def face_keys(self, face: str):
        return self.deck.keys_with_prefix(face + ".") + self.deck.keys_with_prefix("erf.%s." % face)

    def run(self) -> None:
        c, d = self.ctx, self.deck
        if c.periodic[2]:
            self.warn("geometry.is_periodic", "the domain is periodic in z",
                      "A vertical direction with gravity is almost never periodic; set zlo/zhi types.")
        for idim, dim in enumerate("xyz"):
            for side in ("lo", "hi"):
                face = "%s%s" % (dim, side)
                t = (self.fstr(face, "type", "") or "").lower().strip('"')
                keys = self.face_keys(face)
                if c.periodic[idim]:
                    if keys:
                        self.warn(face + ".type", "%s is periodic, so every %s.* key is ignored" % (dim, face),
                                  "Remove the %s.* keys or make the direction non-periodic." % face)
                    continue
                if not t:
                    if c.use_real_bcs and idim < 2:
                        continue
                    self.error(face + ".type", "direction %s is not periodic and %s.type is not set; ERF stops" % (dim, face))
                    continue
                if t not in FACE_TYPES:
                    self.error(face + ".type", "%r is not a boundary type" % t,
                               "Choose from: " + ", ".join(sorted(FACE_TYPES - {"most"})) + ".")
                    continue
                self._face(face, t, idim)
        self._surface_layer()

    def _face(self, face, t, idim) -> None:
        d = self.deck
        if t == "inflow":
            has = any(self.fkey(face, k) for k in ("velocity", "dirichlet_file", "inflow_profile"))
            if not has and not self.bool_("erf.input_bndry_planes", False):
                self.error(face + ".type", "inflow face without %s.velocity, dirichlet_file or inflow_profile; ERF stops" % face)
            vk = self.fkey(face, "velocity")
            v = self.reals(vk) if vk else []
            if v and len(v) != 3:
                self.error(vk, "needs three values")
            if self.fkey(face, "dirichlet_file") and self.fkey(face, "inflow_profile"):
                self.error(self.fkey(face, "dirichlet_file"), "cannot be combined with %s.inflow_profile" % face)
        if t in ("surface_layer", "most") and face != "zlo":
            self.warn(face + ".type", "surface_layer is meant for zlo")
        if t == "surface_layer" and face == "zlo" and self.ctx.is_anelastic is False and idim == 2:
            pass
        for k in ("density", "theta", "velocity"):
            kk = self.fkey(face, k)
            if kk and t in ("outflow", "open", "ho_outflow", "symmetry"):
                self.warn(kk, "not read on a %s face" % t, "Remove the key.")
        if self.fkey(face, "density") and self.fkey(face, "density_grad"):
            self.error(self.fkey(face, "density"), "may not be combined with %s.density_grad" % face)
        if self.fkey(face, "theta") and self.fkey(face, "theta_grad"):
            self.error(self.fkey(face, "theta"), "may not be combined with %s.theta_grad" % face)
        if self.fkey(face, "density") and self.ctx.is_anelastic:
            self.error(self.fkey(face, "density"), "rejected for anelastic runs")
        if t == "outflow" and idim < 2 and not self.ctx.use_real_bcs:
            self.info(face + ".type", "plain outflow reflects acoustic waves; open or inflow_outflow faces are usually quieter")

    def _surface_layer(self) -> None:
        c, d = self.ctx, self.deck
        zlo = (self.fstr("zlo", "type", "") or "").lower()
        has_sl = zlo in ("surface_layer", "most")
        most_keys = d.keys_with_prefix("erf.most.") + d.keys_with_prefix("erf.surface_layer.")
        # per-face surface layers (cloud chamber walls): erf.<face>.most.* and <face>.z0_m etc.
        per_face = [k for k in d.keys() if ".most." in k and not k.startswith("erf.most.")] + \
                   [k for k in d.keys() if k.split(".")[-1] in ("z0_m", "z0_h", "z0_q", "wall_transfer_model")]
        pbl = [p.lower() for p in c.per_level_str("erf.pbl_type", "None")]
        if not has_sl:
            if most_keys:
                self.warn(most_keys[0], "%d erf.most.* / erf.surface_layer.* keys are set but zlo.type is not surface_layer, so they are never read" % len(most_keys),
                          "Set zlo.type = surface_layer, or remove the keys.")
            if any(p in COLUMN_PBL for p in pbl):
                self.error("erf.pbl_type", "the PBL scheme needs zlo.type = surface_layer; ERF stops")
            return
        flux = (self.str_("erf.surface_layer.flux_type", "moeng") or "moeng").lower()
        if flux == "bulk_coeff":
            for k in ("erf.most.Cd",):
                if not d.has(k):
                    self.error(k, "required with flux_type = bulk_coeff")
            if not d.has("erf.most.Ch") and not d.has("erf.most.surf_temp_flux"):
                self.error("erf.most.Ch", "required with flux_type = bulk_coeff unless erf.most.surf_temp_flux is given")
        elif flux in ("custom", "rico"):
            for k in ("erf.most.ustar", "erf.most.tstar"):
                if not d.has(k) and flux == "custom":
                    self.error(k, "required with flux_type = custom")
        if (d.has("erf.most.surf_temp") or d.has("erf.most.pot_surf_temp")) and d.has("erf.most.surf_temp_flux") and flux == "moeng":
            self.warn("erf.most.surf_temp_flux", "both a surface temperature and a heat flux are given; with the Moeng flux type the flux is used and surf_temp only sets the initial surface value",
                      "Keep the one you mean.")
        if per_face:
            self.info(per_face[0], "per-face surface-layer keys are set; the roughness check uses erf.most.z0 only")
            return
        z0 = self.real("erf.most.z0", 0.1)
        if z0 is not None and z0 <= 0 and not d.has("erf.most.roughness_file_name"):
            self.error("erf.most.z0", "must be > 0")
        if z0 and c.dz_min:
            zref = self.real("erf.most.zref", 0.5 * c.dz_min)
            if zref is not None and zref <= z0:
                self.error("erf.most.z0", "roughness z0 = %g m is not below the MOST reference height %g m (log law undefined)" % (z0, zref),
                           "Lower z0 or use a thinner first cell / larger erf.most.zref.")
            elif zref and zref < 5 * z0:
                self.advise("erf.most.z0", "z0 = %g m is within a factor 5 of the reference height %g m; the log law is poorly resolved" % (z0, zref),
                            "Most ABL setups keep the first cell centre at least 10 z0 above the surface.")
        if (self.str_("erf.molec_diff_type", "None") or "None").lower() == "none" and not self._eddy():
            self.warn("zlo.type", "surface_layer with no diffusion operator (molec_diff_type, les_type, rans_type and pbl_type all None): the surface fluxes enter through the diffusive stress, so they have no effect",
                      "Add a closure, or erf.molec_diff_type = ConstantAlpha with zero coefficients if only the surface fluxes are wanted.")

    def _eddy(self) -> bool:
        for key in ("erf.les_type", "erf.rans_type", "erf.pbl_type"):
            if any(v.lower() != "none" for v in self.ctx.per_level_str(key, "None")):
                return True
        return False
