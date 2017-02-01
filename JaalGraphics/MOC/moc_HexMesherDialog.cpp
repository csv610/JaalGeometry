/****************************************************************************
** Meta object code from reading C++ file 'HexMesherDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.5.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "HexMesherDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'HexMesherDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.5.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JHexMesherDialog_t {
    QByteArrayData data[13];
    char stringdata0[231];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JHexMesherDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JHexMesherDialog_t qt_meta_stringdata_JHexMesherDialog = {
    {
        QT_MOC_LITERAL(0, 0, 16), // "JHexMesherDialog"
        QT_MOC_LITERAL(1, 17, 11), // "closeDialog"
        QT_MOC_LITERAL(2, 29, 0), // ""
        QT_MOC_LITERAL(3, 30, 16), // "getGeodeTemplate"
        QT_MOC_LITERAL(4, 47, 19), // "openPolyCubesDialog"
        QT_MOC_LITERAL(5, 67, 18), // "openLegoMeshDialog"
        QT_MOC_LITERAL(6, 86, 20), // "openBernHexOpsDialog"
        QT_MOC_LITERAL(7, 107, 19), // "openMeshStackDialog"
        QT_MOC_LITERAL(8, 127, 20), // "allTet2HexConversion"
        QT_MOC_LITERAL(9, 148, 24), // "openStructuredMeshDialog"
        QT_MOC_LITERAL(10, 173, 20), // "getSchneiderTemplate"
        QT_MOC_LITERAL(11, 194, 13), // "getOctreeMesh"
        QT_MOC_LITERAL(12, 208, 22) // "openSphHexMesherDialog"

    },
    "JHexMesherDialog\0closeDialog\0\0"
    "getGeodeTemplate\0openPolyCubesDialog\0"
    "openLegoMeshDialog\0openBernHexOpsDialog\0"
    "openMeshStackDialog\0allTet2HexConversion\0"
    "openStructuredMeshDialog\0getSchneiderTemplate\0"
    "getOctreeMesh\0openSphHexMesherDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JHexMesherDialog[] = {

// content:
    7,       // revision
    0,       // classname
    0,    0, // classinfo
    11,   14, // methods
    0,    0, // properties
    0,    0, // enums/sets
    0,    0, // constructors
    0,       // flags
    0,       // signalCount

// slots: name, argc, parameters, tag, flags
    1,    0,   69,    2, 0x08 /* Private */,
    3,    0,   70,    2, 0x08 /* Private */,
    4,    0,   71,    2, 0x08 /* Private */,
    5,    0,   72,    2, 0x08 /* Private */,
    6,    0,   73,    2, 0x08 /* Private */,
    7,    0,   74,    2, 0x08 /* Private */,
    8,    0,   75,    2, 0x08 /* Private */,
    9,    0,   76,    2, 0x08 /* Private */,
    10,    0,   77,    2, 0x08 /* Private */,
    11,    0,   78,    2, 0x08 /* Private */,
    12,    0,   79,    2, 0x08 /* Private */,

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

    0        // eod
};

void JHexMesherDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JHexMesherDialog *_t = static_cast<JHexMesherDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0:
            _t->closeDialog();
            break;
        case 1:
            _t->getGeodeTemplate();
            break;
        case 2:
            _t->openPolyCubesDialog();
            break;
        case 3:
            _t->openLegoMeshDialog();
            break;
        case 4:
            _t->openBernHexOpsDialog();
            break;
        case 5:
            _t->openMeshStackDialog();
            break;
        case 6:
            _t->allTet2HexConversion();
            break;
        case 7:
            _t->openStructuredMeshDialog();
            break;
        case 8:
            _t->getSchneiderTemplate();
            break;
        case 9:
            _t->getOctreeMesh();
            break;
        case 10:
            _t->openSphHexMesherDialog();
            break;
        default:
            ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JHexMesherDialog::staticMetaObject = {
    {   &QDialog::staticMetaObject, qt_meta_stringdata_JHexMesherDialog.data,
        qt_meta_data_JHexMesherDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR
    }
};


const QMetaObject *JHexMesherDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JHexMesherDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JHexMesherDialog.stringdata0))
        return static_cast<void*>(const_cast< JHexMesherDialog*>(this));
    if (!strcmp(_clname, "Ui::HexMesherDialog"))
        return static_cast< Ui::HexMesherDialog*>(const_cast< JHexMesherDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JHexMesherDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 11)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 11;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 11)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 11;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
