from __future__ import division, print_function, absolute_import, unicode_literals
#*****************************************************************
#    pyGSTi 0.9:  Copyright 2015 Sandia Corporation
#    This Software is released under the GPL license detailed
#    in the file "license.txt" in the top-level pyGSTi directory
#*****************************************************************
""" Automatic report title generation. """

import numpy as _np

def generate_name():
    adj = _adjectives[ _np.random.randint(0,len(_adjectives)) ]
    noun = _nouns[ _np.random.randint(0,len(_nouns)) ]
    return adj + " " + noun

_nouns = [
    "goatskin",
    "wealth",
    "troll",
    "serpent",
    "sapphire",
    "cleaver",
    "chariot",
    "wizard",
    "pilgrim",
    "carrot",
    "husband",
    "catastrophe",
    "complexity",
    "barbarian",
    "incantation",
    "combat",
    "stewardship",
    "soup",
    "empire",
    "parasite",
    "darkling",
    "sloth",
    "bears",
    "menfolk",
    "dragon",
    "sausage",
    "splinter",
    "honeysuckle",
    "moonbeams",
    "beast",
    "informer",
    "bodyguard",
    "trappings",
    "humanitarian",
    "knife",
    "throne",
    "shaman",
    "hostility",
    "dignity",
    "drapery",
    "spittle",
    "priestess",
    "walnut",
    "plague",
    "twig",
    "devotion",
    "scribe",
    "humans",
    "rival",
    "weapon",
    "curses",
    "hierarchy",
    "enemies",
    "nobleman",
    "weapon",
    "mistress",
    "nostril",
    "fortune",
    "obscenity",
    "partridge",
    "tapestries",
    "putrescence",
    "gossipmonger",
    "orb",
    "haze",
    "provision",
    "shackle" ] + \
    [
    "factory reset button",
    "blood rage",
    "idiot",
    "toaster",
    "legend",
    "death wish",
    "therapy",
    "goal in life",
    "marketing idea",
    "psychic",
    "knife",
    "sandwich",
    "hunting ground",
    "lettuce",
    "kitty",
    "friendly grandma",
    "french chef",
    "corn cake",
    "candlestick maker",
    "coffee pot",
    "brethren",
    "tank",
    "useless brakes",
    "sound barrier",
    "private investor",
    "main people",
    "stock car",
    "elastic band",
    "telephone",
    "mad cow disease",
    "rough-skinned newt",
    "karate",
    "pistol",
    "legal warrant",
    "place of business",
    "double fault",
    "kitty cat",
    "famous landscape painting",
    "hairy legs",
    "old irish cottage",
    "pocket flask",
    "liquid oxygen",
    "laser beams",
    "preventive strike",
    "dingle berry",
    "reading party",
    "generalized bond",
    "indirect expression",
    "messiness",
    "trust fund",
    "volcanic crater",
    "travel guidebook",
    "electric furnace",
    "internal respiration",
    "police squad",
    "mad-dog skullcap",
    "sneaky criminal",
    "keepsake machete",
    "gaming laptop",
    "hissy fit",
    "dog poop",
    "dragon",
    "mediation",
    "patrolman",
    "toilet seat",
    "haunted graveyard",
    "really tough guy",
    "twinkling uncleanness",
    "wrinkle",
    "personal credit line",
    "multi-billionaire",
    "mental disorder",
    "sweet tailpipe",
    "boiling water",
    "deer antler",
    "background story",
    "mood",
    "nibblets",
    "striped hyena",
    "weed whacker" ]

_adjectives = [
    "dead",
    "hairless",
    "sadistic",
    "metal",
    "wild",
    "domesticated",
    "abnormal",
    "medicated",
    "cocky",
    "massive",
    "disrespectful",
    "impressive",
    "out of control",
    "internet worthy",
    "hilarious",
    "tactful",
    "bearded",
    "duck-like",
    "violent",
    "slimy",
    "insanely creepy",
    "talking",
    "angry",
    "shaky",
    "deep",
    "sick",
    "zippy",
    "sticky",
    "fluffy",
    "frozen",
    "filthy",
    "fighting",
    "bonkers",
    "harsh",
    "frisky",
    "greedy",
    "crawly",
    "insane",
    "hideous",
    "ungodly",
    "abusive",
    "drunken",
    "hateful",
    "idiotic",
    "twisted",
    "useless",
    "yapping",
    "magical",
    "confused",
    "flirting",
    "high-end",
    "insecure",
    "maniacal",
    "sickened",
    "slippery",
    "stubborn",
    "vengeful",
    "sinister",
    "costumed",
    "cowardly",
    "haunting",
    "startled",
    "demanding",
    "shivering",
    "offensive",
    "nighttime",
    "startling",
    "disgusting",
    "slap happy",
    "disturbing",
    "blathering",
    "flickering",
    "rebellious",
    "impertinent",
    "bull headed",
    "hyperactive",
    "infuriating",
    "outnumbered",
    "pea-brained",
    "territorial",
    "underhanded",
    "zombie-like",
    "mischievous",
    "at-the-ready",
    "free-loading",
    "house-broken",
    "misunderstood" ]
