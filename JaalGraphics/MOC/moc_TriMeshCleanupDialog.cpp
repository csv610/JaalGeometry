/****************************************************************************
** Meta object code from reading C++ file 'TriMeshCleanupDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "TriMeshCleanupDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'TriMeshCleanupDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JTrimeshCleanupDialog_t {
    QByteArrayData data[17];
    char stringdata0[277];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JTrimeshCleanupDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JTrimeshCleanupDialog_t qt_meta_stringdata_JTrimeshCleanupDialog = {
    {
QT_MOC_LITERAL(0, 0, 21), // "JTrimeshCleanupDialog"
QT_MOC_LITERAL(1, 22, 12), // "makeDelaunay"
QT_MOC_LITERAL(2, 35, 0), // ""
QT_MOC_LITERAL(3, 36, 15), // "displayDelaunay"
QT_MOC_LITERAL(4, 52, 23), // "displayNonManifoldEdges"
QT_MOC_LITERAL(5, 76, 21), // "displayIrregularNodes"
QT_MOC_LITERAL(6, 98, 20), // "displayCircumCircles"
QT_MOC_LITERAL(7, 119, 19), // "displayCollapsables"
QT_MOC_LITERAL(8, 139, 17), // "displayFlippables"
QT_MOC_LITERAL(9, 157, 15), // "openLloydDialog"
QT_MOC_LITERAL(10, 173, 9), // "flipEdges"
QT_MOC_LITERAL(11, 183, 14), // "subdivideEdges"
QT_MOC_LITERAL(12, 198, 11), // "removeEdges"
QT_MOC_LITERAL(13, 210, 17), // "below5DegreeNodes"
QT_MOC_LITERAL(14, 228, 17), // "above7DegreeNodes"
QT_MOC_LITERAL(15, 246, 18), // "openAdvancingfront"
QT_MOC_LITERAL(16, 265, 11) // "closeDialog"

    },
    "JTrimeshCleanupDialog\0makeDelaunay\0\0"
    "displayDelaunay\0displayNonManifoldEdges\0"
    "displayIrregularNodes\0displayCircumCircles\0"
    "displayCollapsables\0displayFlippables\0"
    "openLloydDialog\0flipEdges\0subdivideEdges\0"
    "removeEdges\0below5DegreeNodes\0"
    "above7DegreeNodes\0openAdvancingfront\0"
    "closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JTrimeshCleanupDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      15,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   89,    2, 0x08 /* Private */,
       3,    0,   90,    2, 0x08 /* Private */,
       4,    0,   91,    2, 0x08 /* Private */,
       5,    0,   92,    2, 0x08 /* Private */,
       6,    0,   93,    2, 0x08 /* Private */,
       7,    0,   94,    2, 0x08 /* Private */,
       8,    0,   95,    2, 0x08 /* Private */,
       9,    0,   96,    2, 0x08 /* Private */,
      10,    0,   97,    2, 0x08 /* Private */,
      11,    0,   98,    2, 0x08 /* Private */,
      12,    0,   99,    2, 0x08 /* Private */,
      13,    0,  100,    2, 0x08 /* Private */,
      14,    0,  101,    2, 0x08 /* Private */,
      15,    0,  102,    2, 0x08 /* Private */,
      16,    0,  103,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void JTrimeshCleanupDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JTrimeshCleanupDialog *_t = static_cast<JTrimeshCleanupDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->makeDelaunay(); break;
        case 1: _t->displayDelaunay(); break;
        case 2: _t->displayNonManifoldEdges(); break;
        case 3: _t->displayIrregularNodes(); break;
        case 4: _t->displayCircumCircles(); break;
        case 5: _t->displayCollapsables(); break;
        case 6: _t->displayFlippables(); break;
        case 7: _t->openLloydDialog(); break;
        case 8: _t->flipEdges(); break;
        case 9: _t->subdivideEdges(); break;
        case 10: _t->removeEdges(); break;
        case 11: _t->below5DegreeNodes(); break;
        case 12: _t->above7DegreeNodes(); break;
        case 13: _t->openAdvancingfront(); break;
        case 14: _t->closeDialog(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JTrimeshCleanupDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JTrimeshCleanupDialog.data,
      qt_meta_data_JTrimeshCleanupDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JTrimeshCleanupDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JTrimeshCleanupDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JTrimeshCleanupDialog.stringdata0))
        return static_cast<void*>(const_cast< JTrimeshCleanupDialog*>(this));
    if (!strcmp(_clname, "Ui::TriMeshCleanupDialog"))
        return static_cast< Ui::TriMeshCleanupDialog*>(const_cast< JTrimeshCleanupDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JTrimeshCleanupDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 15)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 15;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 15)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 15;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
